"""
合并的分析流程 (Table A, B, C)
步骤 1 (Table A): 线性回归 (OLS)
    - 暴露 (2): BioX_Adjusted, BioX_Delta
    - 结局: 蛋白质 (2920)
    - 协变量: Age, Sex
    - 输出: 2 个 CSV (TableA_BioX_Adjusted_*.csv, TableA_BioX_Delta_*.csv)
步骤 2 (Table B): Cox 比例风险模型 (CoxPH)
    - 暴露: 蛋白质 (来自 Table A 显著结果, P_bonf < 0.05)
    - 结局 (4): all_cause, cvd, cancer, diabetes
    - 协变量: Age, Sex
    - 输入: 2 个 Table A CSV
    - 输出: 8 个 CSV (TableB_BioX_Adjusted_vs_all_cause_*.csv, ...)
步骤 3 (Table C): 中介分析 (Quasi-Bayesian)
    - 暴露 (2): BioX_Adjusted, BioX_Delta
    - 中介: 蛋白质 (来自 Table B 显著结果, P_bonf < 0.05)
    - 结局 (4): all_cause, cvd, cancer, diabetes (二元事件)
    - 协变量: Age, Sex
    - 输入: 8 个 Table B CSV
    - 输出: 8 个 CSV (TableC_BioX_Adjusted_Mediation_all_cause_*.csv, ...)

总共 18 个输出 CSV 文件。
"""

import os, sys, argparse
import re
import hashlib

# ===================================================================
# 1. 并行配置 (在导入 numpy/scipy/statsmodels 之前)
# ===================================================================

def _preparse_int(argv, keys, default):
    """ 辅助函数：在 argparse 解析前从 sys.argv 提取整数 """
    val = None
    for i, a in enumerate(argv):
        if a in keys and i+1 < len(argv):
            try:
                val = int(argv[i+1])
            except Exception:
                pass
    return val if val is not None else default

# 硬编码并行参数 (不再使用 argparse)
# DEFAULT_JOBS: 使用 90% 的 CPU 核心进行 joblib 并行
DEFAULT_JOBS = max(1, int(0.3 * (os.cpu_count() or 4)))
# DEFAULT_BLAS_THREADS: 每个 joblib 进程只使用 1 个 BLAS 线程，防止超售
DEFAULT_BLAS_THREADS = 1

# 检查用户是否在命令行覆盖了 'jobs' (虽然我们硬编码，但保留此逻辑以防万一)
# 脚本主要依赖硬编码值
JOBS = _preparse_int(sys.argv, ("--jobs","-j"), DEFAULT_JOBS)
BLAS_THREADS = _preparse_int(sys.argv, ("--blas-threads",), DEFAULT_BLAS_THREADS)

# 设置环境变量，强制所有 BLAS (MKL, OpenBLAS) 使用单线程
for var in ("MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[var] = str(BLAS_THREADS)

# 现在可以安全导入
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from lifelines import CoxPHFitter
from joblib import Parallel, delayed
import joblib
from contextlib import contextmanager

from statsmodels.stats.mediation import Mediation
from statsmodels.duration.hazard_regression import PHReg as _SM_PHReg

# 确保 patsy (用于 OLS) 被导入
from patsy import bs  # noqa: F401

# Tqdm (进度条) + Joblib 集成
try:
    from tqdm.auto import tqdm
except ImportError:
    # Fix: Define a dummy class that handles both iterables and manual 'total' usage
    class tqdm:
        def __init__(self, iterable=None, total=None, **kwargs):
            self.iterable = iterable
            self.total = total
        
        def __iter__(self):
            # If used as a loop wrapper: for i in tqdm(range(...))
            if self.iterable is not None:
                yield from self.iterable
            else:
                return iter([])

        def update(self, n=1):
            # If used with joblib: tqdm_object.update()
            pass
        
        def close(self):
            pass
            
        def __enter__(self):
            return self
            
        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

@contextmanager
def tqdm_joblib(tqdm_object):
    """ Tqdm 进度条的 Joblib 回调 """
    class TqdmBatchCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)
    old = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old
        tqdm_object.close()

# ===================================================================
# 2. 全局配置 (硬编码)
# ===================================================================

# ----- 文件与列 -----
INPUT_FILE = "/home/data/laiyuxuan/wc_pro_Batch.csv"
INPUT_SEP = ","
PROTEIN_START_COL_INDEX_1_BASED = 2
PROTEIN_END_COL_INDEX_1_BASED = 2921

# ----- 暴露变量 -----
#WC_VARS = ["BioX_Adjusted", "BioX_Delta"]
WC_VARS = ["WC"]

# ----- 协变量 (与 mediationWC.py 保持一致) -----
# BMI/no_BMI
CONTINUOUS_COVS = ['Age', 'TDI', 'Time_TV', 'MET Physical activity', 'sample_age_days','BMI']
DISCRETE_COVS = [
    'Sex', 'Smoking status', 'Alcohol intake frequency',
    'Oily fish intake', 'Processed meat intake', 'Education',
    'Fruit intake', 'vegetable intake', 'Red meat intake',
    'Ethnic background', 'UK Biobank assessment centre | Instance 0','Batch'
]
ALL_COVS = CONTINUOUS_COVS + DISCRETE_COVS

# ----- 结局/生存定义 -----
DEATH_ICD_COL = "newp_s_alldead"
DEATH_DATE_COL = "new2025516_dead_data"
BASELINE_COL = "date_attending_assessment_centre"
CENSOR_DATE = pd.to_datetime("2024-07-08")
OUTCOME_TYPES = ['all_cause', 'cvd', 'cancer', 'diabetes']
EVENT_COL = "event" # 用于中介分析的二元结局列

# ----- 中介分析参数 (Step 3) -----
MEDIATION_SIMS = 2000
MEDIATION_ALPHA = 0.05
MEDIATION_SEED = 42
MEDIATION_TREAT_MODE = "sd"
MEDIATION_TREAT_DELTA = 1.0
MEDIATION_MIN_EVENTS = 10
MEDIATION_SCALE = True


# ===================================================================
# 3. 辅助函数 (列名清洗、数据处理、ICD 匹配)
# ===================================================================

def clean_col_name(name):
    """ 清洗列名，移除特殊字符 """
    if pd.isna(name): return "Unknown"
    name = str(name).strip()
    name = re.sub(r'[ \(\)\-\/\.,]', '_', name)
    name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')
    return name

def _make_unique(names):
    """ 确保列名唯一 """
    counts = {}
    out = []
    for n in names:
        k = counts.get(n, 0) + 1
        counts[n] = k
        out.append(n if k == 1 else f"{n}_{k}")
    return out

def sanitize_columns_list(cols):
    """ 清洗列名列表并确保唯一 """
    cleaned = [clean_col_name(c) for c in cols]
    return _make_unique(cleaned)

def sanitize_df_columns(df):
    """ 清洗DataFrame的列名 """
    df.columns = sanitize_columns_list(list(df.columns))
    return df

def resolve_name_from_header(raw_name, header_raw, header_sanitized_unique):
    """ 从原始列名解析清洗后的列名 """
    try:
        idx = header_raw.index(raw_name)
        return header_sanitized_unique[idx]
    except Exception:
        return clean_col_name(raw_name)

def apply_full_header_sanitized(df, header_raw, header_sanitized_unique):
    """ 应用完整的列名清洗 """
    df.columns = [resolve_name_from_header(c, header_raw, header_sanitized_unique) for c in df.columns]
    return df

def find_site_column(df):
    """ 查找评估中心列 """
    candidates_raw = ['UK Biobank assessment centre | Instance 0']
    base = clean_col_name(candidates_raw[0])
    if base in df.columns:
        return base
    for c in df.columns:
        if c == base or c.startswith(base + "_"):
            return c
    return None

def _find_by_base(df, base_raw):
    """ 通过基础名称查找列 """
    base = clean_col_name(base_raw)
    if base in df.columns:
        return base
    for c in df.columns:
        if c == base or c.startswith(base + "_"):
            return c
    return None

def impute_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    - Continuous: Median fill
    - TDI: Median fill by Assessment Centre (if available), else global median
    - BMI: Median fill by Sex (if available), else global median
    - Discrete: Mode fill
    - Ethnic background: Fill with 1
    """
    df = df.copy()

    # ---- Continuous: TDI special + others median ----
    targets = [c for c in CONTINUOUS_COVS if c in df.columns]
    for cov in targets:
        tdi_base = clean_col_name('TDI')
        if cov == tdi_base or cov.startswith(tdi_base + "_"):
            df[cov] = pd.to_numeric(df[cov], errors='coerce')
            site_col = find_site_column(df)
            if site_col and site_col in df.columns:
                site_meds = df.groupby(site_col, dropna=False)[cov].median().to_dict()
                overall = df[cov].median()
                def _fill_tdi(row):
                    v = row[cov]
                    if pd.isna(v):
                        return site_meds.get(row.get(site_col), overall)
                    return v
                df[cov] = df.apply(_fill_tdi, axis=1)
            else:
                df[cov] = df[cov].fillna(df[cov].median())
        else:
            series = pd.to_numeric(df[cov], errors='coerce')
            med = series.median()
            if pd.isna(med): med = 0
            df[cov] = series.fillna(med)

    # ---- BMI: median by Sex ----
    bmi_col = _find_by_base(df, 'BMI')
    sex_col = _find_by_base(df, 'Sex')
    if bmi_col and bmi_col in df.columns:
        df[bmi_col] = pd.to_numeric(df[bmi_col], errors='coerce')
        if sex_col and sex_col in df.columns:
            tmp = df[[bmi_col, sex_col]].copy()
            tmp[sex_col] = pd.to_numeric(tmp[sex_col], errors='coerce')
            bmi_medians = tmp.groupby(sex_col, dropna=False)[bmi_col].median()
            overall_median = tmp[bmi_col].median()
            fill_vals = tmp[sex_col].map(bmi_medians).fillna(overall_median)
            df[bmi_col] = df[bmi_col].fillna(fill_vals)
        else:
            df[bmi_col] = df[bmi_col].fillna(df[bmi_col].median())

    # ---- Discrete ----
    eth_col = _find_by_base(df, 'Ethnic background')
    if eth_col and eth_col in df.columns:
        df[eth_col] = df[eth_col].fillna(1)

    for cov in DISCRETE_COVS:
        cov_s = clean_col_name(cov)
        if cov_s not in df.columns:
            continue
        if df[cov_s].isnull().any():
            try:
                m = df[cov_s].mode(dropna=True)
                mode_val = m.iloc[0] if len(m) > 0 else df[cov_s].dropna().iloc[0]
            except Exception:
                mode_val = 0
            df[cov_s] = df[cov_s].fillna(mode_val)

    return df

def parse_first_date(s):
    """ 解析 | 分隔的日期字符串，返回最早的有效日期 """
    if pd.isna(s) or str(s).strip() == "":
        return pd.NaT
    parts = [p.strip() for p in str(s).split("|") if p.strip() != ""]
    if not parts:
        return pd.NaT
    dates = pd.to_datetime(parts, errors="coerce", dayfirst=False, yearfirst=True, infer_datetime_format=True)
    dates = dates.dropna()
    return dates.min() if len(dates) else pd.NaT

def is_event_match(code_str: str, outcome_type: str) -> bool:
    """ 
    检查单个 ICD-10 编码是否匹配指定的结局类型。
    新逻辑: 匹配基于前三位字符 (e.g., I251 视为 I25)。
    """
    if not code_str:
        return False
    
    code = code_str.strip().upper().replace(".", "") # C34.9 -> C349
    if not code: # 处理 " " 或 "."
        return False
        
    # 获取前三位字符作为匹配前缀
    prefix = code[:3] # 'I251' -> 'I25', 'C18' -> 'C18', 'G45' -> 'G45'
    if not prefix:
        return False

    try:
        if outcome_type == 'all_cause':
            # 此逻辑现在在 build_survival_event 中处理
            # 为安全起见保留，但理论上不应为 all_cause 调用此函数
            return True 

        elif outcome_type == 'cvd':
            # G45
            if prefix == 'G45': return True
            
            # I (排除 I79, I80-I89, I95, I97-I99)
            if prefix.startswith('I'):
                # 检查排除项 (基于3位前缀)
                if prefix == 'I79': return False
                if prefix.startswith('I8'): return False # 排除 I80-I89
                if prefix in ('I95', 'I97', 'I98', 'I99'): return False
                # 属于 'I' 且未被排除
                return True
            return False

        elif outcome_type == 'cancer':
            # C00-D48
            if prefix.startswith('C'):
                return True # 涵盖 C00-C97 (e.g., C18, C79)
            
            if prefix.startswith('D'):
                # D00-D48
                # 鲁棒地解析数字 (e.g., D1, D37)
                p_num_str = ""
                for char in prefix[1:]:
                    if char.isdigit():
                        p_num_str += char
                    else:
                        break # 遇到非数字 (e.g., 'DXX')
                
                if not p_num_str: return False # 像 'D' or 'DXX'
                
                p_num = int(p_num_str)
                if len(p_num_str) == 1: # D0-D9
                    p_num = p_num 
                else: # D10-D48
                    p_num = int(p_num_str[:2]) # D37 -> 37, D48 -> 48
                
                if 0 <= p_num <= 48:
                    return True
            return False

        elif outcome_type == 'diabetes':
            # E10-E14
            if prefix.startswith('E1'):
                if len(prefix) < 3: # 'E1'
                    return False # 必须是 E10-E14
                p_num = int(prefix[1:3]) # E10 -> 10, E14 -> 14
                if 10 <= p_num <= 14:
                    return True
            return False
            
    except Exception:
        # 捕获解析 ICD 时的错误 (e.g., 'DXX')
        return False
    
    return False

def build_survival_event(df_in: pd.DataFrame, outcome_type: str) -> pd.DataFrame:
    """
    (用于 Step 2 Cox)
    构建生存数据 (time, event)，正确匹配因果特异性死亡。
    """
    df = df_in.copy()
    base_dt = pd.to_datetime(df[BASELINE_COL], errors="coerce", dayfirst=False, yearfirst=True, infer_datetime_format=True)
    
    # --- 修复: 'all_cause' 使用原始的、能工作的逻辑 ---
    if outcome_type == 'all_cause':
        death_dt = df[DEATH_DATE_COL].apply(parse_first_date)
        has_icd = (~df[DEATH_ICD_COL].astype(str).fillna("").eq(""))
        died = has_icd & death_dt.notna() & (death_dt <= CENSOR_DATE)
        
        observed_dt = death_dt.copy()
        observed_dt.loc[~died] = CENSOR_DATE
        
    # --- 保持特定原因的逻辑 ---
    else:
        icd_codes_series = df[DEATH_ICD_COL].astype(str).fillna("").str.split('|')
        death_dates_series = df[DEATH_DATE_COL].astype(str).fillna("").str.split('|')

        event_dates = pd.Series(pd.NaT, index=df.index)
        
        for idx in df.index:
            codes = icd_codes_series.loc[idx]
            dates_str = death_dates_series.loc[idx]
            
            if not isinstance(codes, list) or not isinstance(dates_str, list):
                continue

            dates = pd.to_datetime(dates_str, errors='coerce', dayfirst=False, yearfirst=True, infer_datetime_format=True)
            
            valid_event_dts = []
            for code, dt in zip(codes, dates):
                # 修复: 使用 .strip() 确保 ' ' 或 '||' 中的空条目被视为空
                code_clean = str(code).strip()
                if pd.isna(dt) or dt > CENSOR_DATE or not code_clean:
                    continue
                
                if is_event_match(code_clean, outcome_type):
                    valid_event_dts.append(dt)
            
            if valid_event_dts:
                event_dates.loc[idx] = min(valid_event_dts)

        # 计算 time 和 event
        died = event_dates.notna()
        observed_dt = event_dates.copy()
        observed_dt.loc[~died] = CENSOR_DATE
    
    # --- 通用逻辑 ---
    diff = (observed_dt - base_dt)
    time_days = diff.dt.total_seconds() / (24.0 * 3600.0)

    out = df.copy()
    out["time"] = time_days
    out["event"] = died.astype(int)
    
    # 过滤无效时间
    out = out[out["time"].notna() & (out["time"] > 0)]
    return out

def build_binary_event(df_in: pd.DataFrame, outcome_type: str) -> pd.DataFrame:
    """
    (用于 Step 3 中介分析)
    构建二元结局 (event=0/1)，匹配因果特异性死亡。
    """
    df = df_in.copy()
    
    # --- 修复: 'all_cause' 使用原始的、能工作的逻辑 ---
    if outcome_type == 'all_cause':
        death_dt = df[DEATH_DATE_COL].apply(parse_first_date)
        has_icd = (~df[DEATH_ICD_COL].astype(str).fillna("").eq(""))
        died = has_icd & death_dt.notna() & (death_dt <= CENSOR_DATE)
        event = died.astype(int)
        
    # --- 保持特定原因的逻辑 ---
    else:
        icd_codes_series = df[DEATH_ICD_COL].astype(str).fillna("").str.split('|')
        death_dates_series = df[DEATH_DATE_COL].astype(str).fillna("").str.split('|')

        event = pd.Series(0, index=df.index, dtype=int)
        
        for idx in df.index:
            codes = icd_codes_series.loc[idx]
            dates_str = death_dates_series.loc[idx]
            
            if not isinstance(codes, list) or not isinstance(dates_str, list):
                continue

            dates = pd.to_datetime(dates_str, errors='coerce', dayfirst=False, yearfirst=True, infer_datetime_format=True)
            
            for code, dt in zip(codes, dates):
                # 修复: 使用 .strip() 确保 ' ' 或 '||' 中的空条目被视为空
                code_clean = str(code).strip()
                if pd.isna(dt) or dt > CENSOR_DATE or not code_clean:
                    continue
                
                if is_event_match(code_clean, outcome_type):
                    event.loc[idx] = 1
                    break # 只要发生过一次就算

    # --- 通用逻辑 ---
    out = df.copy()
    out[EVENT_COL] = event
    return out

# ===================================================================
# 4. Step 1: OLS (WC vs Proteins) - 使用矩阵运算优化
# ===================================================================

def run_step1_ols(input_file, sep, wc_col, protein_indices,
                  all_covs_list, discrete_covs_list, out_file):
    """
    (Step 1 主函数)
    执行 OLS 分析 (WC vs 所有蛋白质)。
    使用矩阵运算优化，与 mediationWC.py 保持一致。
    """
    print(f"  [Step 1] 读取数据: {input_file} (用于 {wc_col})")
    df_raw = pd.read_csv(input_file, sep=sep, dtype=str)
    df_raw = sanitize_df_columns(df_raw)

    wc_col = clean_col_name(wc_col)

    # 提取蛋白列
    pro_start = protein_indices[0] - 1
    pro_end   = protein_indices[1]
    protein_cols = list(df_raw.columns[pro_start:pro_end])
    p_total = len(protein_cols)
    print(f"  [Step 1] 找到 {p_total} 个蛋白质 (索引 {protein_indices[0]}~{protein_indices[1]})")

    # keep covs + site(for TDI imputation)
    keep_cols = [wc_col] + [clean_col_name(c) for c in all_covs_list if clean_col_name(c) in df_raw.columns]
    site_col = find_site_column(df_raw)
    if site_col and site_col in df_raw.columns:
        keep_cols.append(site_col)
    keep_cols = list(dict.fromkeys([c for c in keep_cols if c in df_raw.columns]))

    df_cov = df_raw[keep_cols].copy()

    # --- WC + Covariates to numeric & site kept raw ---
    cont_cols = [c for c in CONTINUOUS_COVS if c in df_cov.columns]
    num_cols = list(dict.fromkeys([wc_col] + cont_cols))
    for c in num_cols:
        df_cov[c] = pd.to_numeric(df_cov[c], errors="coerce")

    # 剔除 WC 缺失
    n_before = len(df_cov)
    df_cov = df_cov[~df_cov[wc_col].isna()].copy()
    print(f"  [Step 1] 剔除 {wc_col} 缺失: {n_before - len(df_cov)} 行, 保留 {len(df_cov)} 行")

    # 协变量插补
    df_cov = impute_data(df_cov)

    # build design matrix for covariates (excluding WC)
    disc_cols = [clean_col_name(c) for c in discrete_covs_list]
    disc_cols = [c for c in disc_cols if c in df_cov.columns]
    cont_cols = [c for c in cont_cols if c in df_cov.columns]

    X_parts = []
    if cont_cols:
        X_parts.append(df_cov[cont_cols].apply(pd.to_numeric, errors="coerce"))
    if disc_cols:
        X_parts.append(pd.get_dummies(df_cov[disc_cols], prefix=disc_cols, drop_first=True))

    if X_parts:
        X_cov = pd.concat(X_parts, axis=1)
    else:
        X_cov = pd.DataFrame(index=df_cov.index)

    X_cov = X_cov.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    n = df_cov.shape[0]
    C = np.column_stack([np.ones(n, dtype=np.float64), X_cov.to_numpy(np.float64)])
    wc = df_cov[wc_col].to_numpy(np.float64)

    try:
        Q, R = np.linalg.qr(C, mode="reduced")
        wc_coef = np.linalg.solve(R, Q.T @ wc)
        wc_hat = C @ wc_coef
    except Exception:
        wc_coef, *_ = np.linalg.lstsq(C, wc, rcond=None)
        wc_hat = C @ wc_coef

    wc_res = wc - wc_hat
    denom = float(np.dot(wc_res, wc_res))
    if not np.isfinite(denom) or denom <= 0:
        print("  [Step 1][ERROR] WC has zero variance after covariate adjustment.")
        out = pd.DataFrame({"protein": protein_cols, "Beta": np.nan, "SE": np.nan, "P_value": np.nan})
        out["P_bonferroni"] = np.nan
        out["Protein_count"] = p_total
        out.to_csv(out_file, index=False)
        return out_file

    p_full = C.shape[1] + 1
    df_resid = n - p_full
    if df_resid <= 1:
        print("  [Step 1][ERROR] Not enough degrees of freedom.")
        df_resid = max(2, df_resid)

    print("  [Step 1] 处理蛋白质...")
    df_pro = df_raw.loc[df_cov.index, protein_cols].apply(pd.to_numeric, errors="coerce")
    medians = df_pro.median(axis=0, skipna=True)
    df_pro = df_pro.fillna(medians).fillna(0.0)

    results = []
    block = 256

    for i in tqdm(range(0, p_total, block), desc=f"Step 1 OLS-matrix ({wc_col})", unit="block"):
        cols_blk = protein_cols[i:i+block]
        Y = df_pro[cols_blk].to_numpy(np.float64)

        try:
            if "Q" in locals():
                B = np.linalg.solve(R, Q.T @ Y)
                Y_hat = C @ B
            else:
                B, *_ = np.linalg.lstsq(C, Y, rcond=None)
                Y_hat = C @ B
        except Exception:
            B, *_ = np.linalg.lstsq(C, Y, rcond=None)
            Y_hat = C @ B

        Y_res = Y - Y_hat

        num = (wc_res[:, None] * Y_res).sum(axis=0)
        beta = num / denom

        resid = Y_res - wc_res[:, None] * beta
        rss = (resid ** 2).sum(axis=0)
        sigma2 = rss / df_resid
        se = np.sqrt(sigma2 / denom)
        tval = beta / se
        pval = 2.0 * stats.t.sf(np.abs(tval), df=df_resid)

        for j, prot in enumerate(cols_blk):
            results.append([prot, float(beta[j]), float(se[j]), float(pval[j])])

    out = pd.DataFrame(results, columns=["protein", "Beta", "SE", "P_value"])
    out["P_bonferroni"] = (out["P_value"] * p_total).clip(upper=1.0)
    out["Protein_count"] = p_total
    out = out.sort_values("P_value").reset_index(drop=True)
    out.to_csv(out_file, index=False)
    print(f"  [Step 1] OLS 结果已保存 -> {out_file}")
    return out_file

# ===================================================================
# 5. Step 2: Cox (Protein vs Outcome) - 与 mediationWC.py 保持一致
# ===================================================================

def run_cox_one(prot, df_base, covar_cols):
    """
    (Step 2 并行单元)
    为单个蛋白质运行 CoxPH。
    返回: [prot, hr, ci_l, ci_u, p]
    """
    try:
        cols = ["time", "event"] + covar_cols + [prot]
        dd = df_base[cols].copy()
        
        # 确保所有列都是数值型
        for c in covar_cols + [prot]:
            dd[c] = pd.to_numeric(dd[c], errors='coerce')
        
        dd = dd.dropna()
        
        if dd.shape[0] < 50 or dd["event"].sum() < 10:
            return [prot, np.nan, np.nan, np.nan, np.nan]
            
        cph = CoxPHFitter()
        cph.fit(dd, duration_col="time", event_col="event", show_progress=False)
                
        if prot not in cph.summary.index:
            return [prot, np.nan, np.nan, np.nan, np.nan]
            
        row = cph.summary.loc[prot]
        # lifelines >= 0.27.0
        hr   = float(row.get('exp(coef)', np.exp(row['coef'])))
        ci_l = float(row.get('exp(coef) lower 95%', np.exp(row['coef lower 95%'])))
        ci_u = float(row.get('exp(coef) upper 95%', np.exp(row['coef upper 95%'])))
        p    = float(row['p'])
        
        return [prot, hr, ci_l, ci_u, p]
        
    except Exception:
        # 捕获 lifelines 拟合失败
        return [prot, np.nan, np.nan, np.nan, np.nan]

def run_step2_cox(input_file, sep, table_a_csv, outcome, all_covs_list, out_file):
    """
    (Step 2 主函数)
    执行 CoxPH 分析 (显著蛋白 vs 单个结局)。
    """
    
    # 1. 候选蛋白 (来自 TableA, P_bonf < 0.05)
    try:
        tableA = pd.read_csv(table_a_csv)
        prot_col = 'protein' if 'protein' in tableA.columns else 'Protein'
        p_col = 'P_bonferroni'
        
        cand_prots = (tableA.loc[tableA[p_col] < 0.05, prot_col]
                              .dropna().astype(str).tolist())
        cand_prots = list(dict.fromkeys(cand_prots)) # 去重
    except Exception as e:
        print(f"  [Step 2] 无法读取 TableA: {table_a_csv}。跳过。错误: {e}")
        return None

    if not cand_prots:
        print(f"  [Step 2] {table_a_csv} 中没有 P_bonf < 0.05 的显著蛋白。跳过 {outcome}。")
        return None
        
    print(f"  [Step 2] 找到 {len(cand_prots)} 个候选蛋白 (来自 {table_a_csv})")

    # 2. 读取主数据
    print(f"  [Step 2] 读取数据: {input_file} (用于 {outcome})")
    df_full = pd.read_csv(input_file, sep=sep, dtype=str)
    df_full = sanitize_df_columns(df_full)
    
    # 3. 构建生存数据 (因果特异性)
    print(f"  [Step 2] 构建 {outcome} 生存数据...")
    df_surv = build_survival_event(df_full, outcome)
    
    # 4. 协变量插补
    df_surv = impute_data(df_surv)

    # 5. 蛋白中位数插补 (仅候选蛋白)
    print(f"  [Step 2] 插补 {len(cand_prots)} 个蛋白质的中位数...")
    present_prots = [p for p in cand_prots if p in df_surv.columns]
    for p in present_prots:
        s = pd.to_numeric(df_surv[p], errors='coerce')
        med = np.nanmedian(s)
        if np.isfinite(med):
            df_surv[p] = np.where(np.isnan(s), med, s)
        else:
            df_surv[p] = s.fillna(0)
            print(f"  [WARN] 蛋白 {p} 的中位数无效，用 0 填充 NA。")
        
    if not present_prots:
        print(f"  [Step 2] 候选蛋白在主数据 {input_file} 中均未找到。跳过。")
        return None

    # 6. 准备基础数据
    sanitized_all_covs = [clean_col_name(c) for c in all_covs_list]
    sanitized_discrete = [clean_col_name(c) for c in DISCRETE_COVS]
    
    keep_cols = set(["time","event"]) | set(sanitized_all_covs) | set(present_prots)
    df_base = df_surv.loc[:, [c for c in df_surv.columns if c in keep_cols]].copy()

    # 7. 处理协变量：离散变量生成哑变量
    final_covar_cols = []
    for c in sanitized_all_covs:
        if c not in df_base.columns:
            continue
        if c in sanitized_discrete:
            dummies = pd.get_dummies(df_base[c], prefix=c, drop_first=True)
            for d_col in dummies.columns:
                df_base[d_col] = dummies[d_col]
                final_covar_cols.append(d_col)
        else:
            df_base[c] = pd.to_numeric(df_base[c], errors='coerce')
            final_covar_cols.append(c)
        
    print(f"  [Step 2] 开始并行 CoxPH (n_jobs={JOBS}) ...")

    # 8. 并行执行
    with tqdm_joblib(tqdm(total=len(present_prots), desc=f"Step 2 Cox ({outcome})", unit="prot")) as progress_bar:
        results = Parallel(n_jobs=JOBS, backend="loky")(
            delayed(run_cox_one)(p, df_base, final_covar_cols) for p in present_prots
        )

    print("  [Step 2] 整理结果...")
    out = pd.DataFrame(results, columns=["protein", "HR", "CI_low", "CI_high", "P_value"])
    
    # Bonferroni (按实际测试的蛋白数 m)
    m = max(1, len(present_prots))
    out["P_bonf"] = (out["P_value"] * m).clip(upper=1.0)
    out["HR_95CI"] = out.apply(
        lambda r: f"{r.HR:.3f} ({r.CI_low:.3f}-{r.CI_high:.3f})" if pd.notna(r.HR) else "",
        axis=1
    )
    
    out = out.sort_values("P_value", na_position="last").reset_index(drop=True)
    out.to_csv(out_file, index=False)
    print(f"  [Step 2] Cox 结果已保存 -> {out_file}")
    
    return out_file

# ===================================================================
# 6. Step 3: Mediation (WC -> Protein -> Outcome) - 使用 statsmodels.Mediation
# ===================================================================

class PHRegPredOnly(_SM_PHReg):
    """
    Override to make PHReg.predict(params, exog) return ndarray directly.
    Statsmodels Mediation.fit expects this, but PHReg usually returns a bunch.
    """
    def predict(self, params, exog=None, cov_params=None, endog=None,
                strata=None, offset=None, pred_type="lhr", pred_only=False):
        return super().predict(
            params, exog=exog, cov_params=cov_params, endog=endog,
            strata=strata, offset=offset, pred_type="lhr", pred_only=True
        )

def _pvalue(vec):
    """ 计算 p 值 """
    vec = np.asarray(vec)
    vec = vec[np.isfinite(vec)]
    if vec.size == 0:
        return np.nan
    return 2 * min(np.sum(vec > 0), np.sum(vec < 0)) / float(len(vec))

def mediate_one_cox_sm(prot, dfb, wc_col, cov_cols_numeric,
                       sims, alpha, seed, treat_mode, treat_delta,
                       min_events, scale_mediator):
    """
    Step3: Cox PHReg + statsmodels Mediation
    Effect scale: log-hazard (lhr) differences.
    """
    need_cols = ["time", "event", wc_col] + cov_cols_numeric + [prot]
    dd = dfb.loc[:, [c for c in need_cols if c in dfb.columns]].copy()
    dd = dd.dropna()

    n_eff = int(dd.shape[0])
    if n_eff < 200 or dd["event"].sum() < min_events:
        return [prot] + [np.nan]*12 + [n_eff]

    # ---- mediator z-score (optional) ----
    M = pd.to_numeric(dd[prot], errors="coerce")
    if scale_mediator:
        mu = M.mean()
        sd = M.std(ddof=1)
        if not np.isfinite(sd) or sd == 0:
            return [prot] + [np.nan]*12 + [n_eff]
        M = (M - mu) / sd

    # ---- exposure scaling so that 0 vs 1 == (mean vs mean+1SD) or (mean vs mean+delta) ----
    X_raw = pd.to_numeric(dd[wc_col], errors="coerce")
    mu_x = X_raw.mean()
    if treat_mode == "delta":
        denom = float(treat_delta)
        if not np.isfinite(denom) or denom == 0:
            return [prot] + [np.nan]*12 + [n_eff]
        X = (X_raw - mu_x) / denom
    else:
        sd_x = X_raw.std(ddof=1)
        if not np.isfinite(sd_x) or sd_x == 0:
            return [prot] + [np.nan]*12 + [n_eff]
        X = (X_raw - mu_x) / sd_x

    # numeric covs
    COVS = dd[cov_cols_numeric].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # time/event
    T = pd.to_numeric(dd["time"], errors="coerce")
    E = pd.to_numeric(dd["event"], errors="coerce")

    # final assemble
    dat = pd.DataFrame({"X": X, "M": M, "time": T, "event": E}, index=dd.index)
    dat = pd.concat([dat, COVS], axis=1).dropna()

    n_eff = int(dat.shape[0])
    if n_eff < 200 or dat["event"].sum() < min_events:
        return [prot] + [np.nan]*12 + [n_eff]

    try:
        # reproducible per-protein
        h = int(hashlib.md5(prot.encode("utf-8")).hexdigest()[:8], 16)
        np.random.seed(seed + (h % 1_000_000))

        # mediator model: OLS with intercept
        mex = sm.add_constant(dat[["X"] + cov_cols_numeric], has_constant="add")
        mediator_model = sm.OLS(dat["M"].to_numpy(np.float64), mex.to_numpy(np.float64))

        # outcome model: PHReg (NO intercept!)   exog order: [X, M, covs]
        oex = dat[["X", "M"] + cov_cols_numeric].to_numpy(np.float64)
        outcome_model = PHRegPredOnly(dat["time"].to_numpy(np.float64), oex,
                                      status=dat["event"].to_numpy(np.int8),
                                      ties="breslow")

        # exposure positions: outcome exog X is col0; mediator exog has const then X so X is col1
        tx_pos = [0, 1]
        # mediator position in outcome exog: M is col1
        med_pos = 1

        med = Mediation(outcome_model, mediator_model, tx_pos, med_pos).fit(
            method="parametric", n_rep=int(sims)
        )

        # Use "average" effects
        ADE_vec  = np.asarray(med.ADE_avg,  dtype=float)
        ACME_vec = np.asarray(med.ACME_avg, dtype=float)
        PM_vec   = np.asarray(med.prop_med_avg, dtype=float)

        # point estimates: ADE/ACME mean; PM median (same as statsmodels summary behavior)
        ADE_hat  = float(np.nanmean(ADE_vec))
        ACME_hat = float(np.nanmean(ACME_vec))
        PM_hat   = float(np.nanmedian(PM_vec))

        def ci(arr, use_median=False):
            arr = np.asarray(arr, dtype=float)
            arr = arr[np.isfinite(arr)]
            if arr.size == 0:
                return np.nan, np.nan, np.nan, np.nan
            lo, hi = np.percentile(arr, [100*alpha/2, 100*(1-alpha/2)])
            p = _pvalue(arr)
            est = float(np.median(arr)) if use_median else float(np.mean(arr))
            return est, float(lo), float(hi), float(p)

        _, ade_lo,  ade_hi,  ade_p  = ci(ADE_vec,  use_median=False)
        _, acme_lo, acme_hi, acme_p = ci(ACME_vec, use_median=False)
        _, pm_lo,   pm_hi,   pm_p   = ci(PM_vec,   use_median=True)

        # total effect = ACME + ADE (on lhr scale, averaged)
        TE_hat = float(np.nanmean(np.asarray(med.total_effect, dtype=float))) if hasattr(med, "total_effect") else np.nan
        if not np.isfinite(TE_hat) or TE_hat == 0:
            # PM already computed from statsmodels ratio; keep it
            pass

        return [prot,
                ADE_hat,  ade_lo,  ade_hi,  ade_p,
                ACME_hat, acme_lo, acme_hi, acme_p,
                PM_hat,   pm_lo,   pm_hi,   pm_p,
                n_eff]
    except Exception:
        return [prot] + [np.nan]*12 + [n_eff]

def run_step3_mediation(input_file, sep, table_b_csv, wc_col, outcome, 
                        all_covs_list, out_file):
    """
    (Step 3 主函数)
    执行中介分析 (WC -> 显著蛋白 -> 单个结局)
    使用 statsmodels.Mediation (与 mediationWC.py 保持一致)
    """
    
    # 1. 候选蛋白 (来自 TableB, P_bonf < 0.05)
    try:
        tableB = pd.read_csv(table_b_csv)
        prot_col = 'protein' if 'protein' in tableB.columns else 'Protein'
        p_col = 'P_bonf'
        
        cand_prots = (tableB.loc[tableB[p_col] < 0.05, prot_col]
                              .dropna().astype(str).tolist())
        cand_prots = list(dict.fromkeys(cand_prots))
    except Exception as e:
        print(f"  [Step 3] 无法读取 TableB: {table_b_csv}。跳过。错误: {e}")
        return None

    if not cand_prots:
        print(f"  [Step 3] {table_b_csv} 中没有 P_bonf < 0.05 的显著蛋白。跳过 {outcome}。")
        return None
        
    print(f"  [Step 3] 找到 {len(cand_prots)} 个候选蛋白 (来自 {table_b_csv})")

    # 2. 读取主数据
    print(f"  [Step 3] 读取数据: {input_file} (用于 {wc_col} -> {outcome})")
    df_full = pd.read_csv(input_file, sep=sep, dtype=str)
    df_full = sanitize_df_columns(df_full)
    wc_col = clean_col_name(wc_col)

    # 3. 构建生存数据 (因果特异性)
    print(f"  [Step 3] 构建 {outcome} 生存数据...")
    df_surv = build_survival_event(df_full, outcome)
    
    # 4. 协变量插补
    df_surv = impute_data(df_surv)

    # 5. 准备基础数据
    present_prots = [p for p in cand_prots if p in df_surv.columns]
    if not present_prots:
        print(f"  [Step 3] 候选蛋白在主数据 {input_file} 中均未找到。跳过。")
        return None

    sanitized_all_covs = [clean_col_name(c) for c in all_covs_list]
    sanitized_discrete = [clean_col_name(c) for c in DISCRETE_COVS]

    # 蛋白质数值化 + 中位数插补
    print(f"  [Step 3] 转换蛋白质为数值并插补...")
    for p in present_prots:
        s = pd.to_numeric(df_surv[p], errors='coerce')
        medv = np.nanmedian(s)
        df_surv[p] = np.where(np.isnan(s), medv, s) if np.isfinite(medv) else s.fillna(0)

    if wc_col not in df_surv.columns:
        print(f"  [Step 3] 暴露变量 {wc_col} 不在数据中。跳过。")
        return None

    keep_cols = ["time", "event", wc_col] + sanitized_all_covs + present_prots
    dfb = df_surv.loc[:, [c for c in keep_cols if c in df_surv.columns]].copy()
    dfb[wc_col] = pd.to_numeric(dfb[wc_col], errors='coerce')

    # 6. 处理协变量：离散变量生成哑变量
    final_covar_cols_numeric = []
    for c in sanitized_all_covs:
        if c not in dfb.columns:
            continue
        if c in sanitized_discrete:
            dummies = pd.get_dummies(dfb[c], prefix=c, drop_first=True)
            for d_col in dummies.columns:
                dfb[d_col] = dummies[d_col]
                final_covar_cols_numeric.append(d_col)
        else:
            dfb[c] = pd.to_numeric(dfb[c], errors='coerce')
            final_covar_cols_numeric.append(c)

    print(f"  [Step 3] 开始并行中介分析 (PHReg + statsmodels.Mediation) (n_jobs={JOBS}) ...")

    # 7. 并行执行
    with tqdm_joblib(tqdm(total=len(present_prots), desc=f"Step 3 Med (Cox) ({wc_col}->{outcome})", unit="prot")) as progress_bar:
        results = Parallel(n_jobs=JOBS, backend="loky")(
            delayed(mediate_one_cox_sm)(
                p, dfb, wc_col, final_covar_cols_numeric,
                sims=MEDIATION_SIMS, alpha=MEDIATION_ALPHA, seed=MEDIATION_SEED,
                treat_mode=MEDIATION_TREAT_MODE, treat_delta=MEDIATION_TREAT_DELTA,
                min_events=MEDIATION_MIN_EVENTS,
                scale_mediator=MEDIATION_SCALE
            )
            for p in present_prots
        )

    print("  [Step 3] 整理结果...")
    cols = ["protein",
            "ADE", "ADE_low", "ADE_high", "ADE_p",
            "ACME","ACME_low","ACME_high","ACME_p",
            "PM",  "PM_low",  "PM_high",  "PM_p",
            "N"]
    out = pd.DataFrame(results, columns=cols)
    
    # Bonferroni 校正 (ACME P 值)
    m = max(1, len(present_prots))
    out["ACME_p_bonf"] = (out["ACME_p"] * m).clip(upper=1.0)
    
    # 格式化
    def fmt_ci(v, lo, hi):
        if not np.isfinite(v) or not np.isfinite(lo) or not np.isfinite(hi): return ""
        return f"{v:.4e} ({lo:.4e},{hi:.4e})"
    def fmt_pm(pm, lo, hi):
        if not np.isfinite(pm) or not np.isfinite(lo) or not np.isfinite(hi): return ""
        return f"{pm*100:.2f}% ({lo*100:.2f}%,{hi*100:.2f}%)"

    out["ADE_str"]  = out.apply(lambda r: fmt_ci(r["ADE"],  r["ADE_low"],  r["ADE_high"]), axis=1)
    out["ACME_str"] = out.apply(lambda r: fmt_ci(r["ACME"], r["ACME_low"], r["ACME_high"]), axis=1)
    out["PM_str"]   = out.apply(lambda r: fmt_pm(r["PM"],   r["PM_low"],   r["PM_high"]), axis=1)

    out = out.sort_values(["ACME_p","ADE_p"], na_position="last").reset_index(drop=True)
    out.to_csv(out_file, index=False)
    print(f"  [Step 3] 中介分析结果已保存 -> {out_file}")
    
    return out_file

# ===================================================================
# 7. 主执行流程 (Main)
# ===================================================================

def main():
    global BASELINE_COL, DEATH_DATE_COL, DEATH_ICD_COL
    global CONTINUOUS_COVS, DISCRETE_COVS, ALL_COVS, WC_VARS

    print("==========================================================")
    print(" 开始执行合并分析流程 (Table A, B, C)")
    print(f" 并行设置: n_jobs={JOBS}, blas_threads={BLAS_THREADS}")
    print(f" 主输入文件: {INPUT_FILE}")
    print(f" 结局类型: {OUTCOME_TYPES}")
    print(f" 协变量: {ALL_COVS}")
    print("==========================================================")

    # 读取原始表头并清洗列名
    header_raw = pd.read_csv(INPUT_FILE, sep=INPUT_SEP, nrows=0).columns.tolist()
    header_sanitized_unique = sanitize_columns_list(header_raw)

    # 解析清洗后的列名
    BASELINE_COL    = resolve_name_from_header(BASELINE_COL,    header_raw, header_sanitized_unique)
    DEATH_DATE_COL  = resolve_name_from_header(DEATH_DATE_COL,  header_raw, header_sanitized_unique)
    DEATH_ICD_COL   = resolve_name_from_header(DEATH_ICD_COL,   header_raw, header_sanitized_unique)

    # 清洗协变量列名
    CONTINUOUS_COVS = [clean_col_name(c) for c in CONTINUOUS_COVS]
    DISCRETE_COVS   = [clean_col_name(c) for c in DISCRETE_COVS]
    ALL_COVS        = CONTINUOUS_COVS + DISCRETE_COVS
    WC_VARS         = [clean_col_name(c) for c in WC_VARS]

    # 新增: 打印结局数量
    print("\n[--- 正在加载数据并计算结局数量 ---]")
    try:
        # 仅加载计算结局所需的列，以节省内存
        cols_to_load_raw = [BASELINE_COL, DEATH_DATE_COL, DEATH_ICD_COL]
        present_raw = [c for c in cols_to_load_raw if c in header_raw]
        df_counts = pd.read_csv(INPUT_FILE, sep=INPUT_SEP, usecols=present_raw, dtype=str)
        df_counts = apply_full_header_sanitized(df_counts, header_raw, header_sanitized_unique)

        print(f"  总加载行数: {len(df_counts)}")
        
        event_counts = {}
        for outcome in OUTCOME_TYPES:
            # 使用 build_binary_event 来计算总事件数
            df_events = build_binary_event(df_counts, outcome)
            count = df_events[EVENT_COL].sum()
            event_counts[outcome] = count
            print(f"  [结局统计] {outcome}: {count} 个事件")
            
    except Exception as e:
        print(f"[ERROR] 无法加载数据或计算结局数量: {e}")
        print(f"详细错误: {str(e)}")
        return # 如果无法加载数据，则停止执行

    # 存储中间文件名
    step1_results = {} # key: wc_var, value: file_path
    step2_results = {} # key: (wc_var, outcome), value: file_path

    protein_indices = (PROTEIN_START_COL_INDEX_1_BASED, PROTEIN_END_COL_INDEX_1_BASED)

    # --- 执行 Step 1 ---
    print("\n[--- 开始 Step 1: OLS (WC vs Proteins) ---]")
    for wc_var in WC_VARS:
        print(f"\n[Step 1] 分析暴露: {wc_var}")
        out_file = f"TableA_{wc_var}_vs_Proteins_OLS.csv"
        try:
            path = run_step1_ols(
                INPUT_FILE, INPUT_SEP, wc_var, protein_indices, 
                ALL_COVS, DISCRETE_COVS, out_file
            )
            step1_results[wc_var] = path
        except Exception as e:
            print(f"[ERROR] Step 1 ({wc_var}) 失败: {e}")
            print(f"详细错误: {str(e)}")

    # --- 执行 Step 2 ---
    print("\n[--- 开始 Step 2: Cox (Protein vs Outcomes) ---]")
    for wc_var in WC_VARS:
        table_a_file = step1_results.get(wc_var)
        if not table_a_file:
            print(f"\n[Step 2] 跳过 {wc_var} (因 Step 1 失败或未运行)")
            continue
        
        print(f"\n[Step 2] 基于 {wc_var} 的结果 (TableA: {table_a_file})")
        
        for outcome in OUTCOME_TYPES:
            print(f"  [Step 2] 分析结局: {outcome}")
            out_file = f"TableB_{wc_var}_vs_{outcome}_Cox.csv"
            try:
                path = run_step2_cox(
                    INPUT_FILE, INPUT_SEP, table_a_file, 
                    outcome, ALL_COVS, out_file
                )
                if path:
                    step2_results[(wc_var, outcome)] = path
            except Exception as e:
                print(f"[ERROR] Step 2 ({wc_var} vs {outcome}) 失败: {e}")
                print(f"详细错误: {str(e)}")

    # --- 执行 Step 3 ---
    print("\n[--- 开始 Step 3: Mediation (WC -> Protein -> Outcome) ---]")
    for wc_var in WC_VARS:
        print(f"\n[Step 3] 分析暴露: {wc_var}")
        
        for outcome in OUTCOME_TYPES:
            table_b_file = step2_results.get((wc_var, outcome))
            if not table_b_file:
                print(f"  [Step 3] 跳过 {wc_var} -> {outcome} (因 Step 2 失败或无显著蛋白)")
                continue

            print(f"  [Step 3] 分析结局: {outcome} (TableB: {table_b_file})")
            out_file = f"TableC_{wc_var}_Mediation_{outcome}.csv"
            try:
                run_step3_mediation(
                    INPUT_FILE, INPUT_SEP, table_b_file,
                    wc_col=wc_var, outcome=outcome, 
                    all_covs_list=ALL_COVS, out_file=out_file
                )
            except Exception as e:
                print(f"[ERROR] Step 3 ({wc_var} -> {outcome}) 失败: {e}")
                print(f"详细错误: {str(e)}")

    print("\n==========================================================")
    print(" 所有分析流程执行完毕。")
    print("==========================================================")

if __name__ == "__main__":
    main()