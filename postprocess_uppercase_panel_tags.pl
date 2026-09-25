#!/usr/bin/perl
## Replace lower-case panel tags with upper-case tags in a one-page figure PDF, keeping the vector content:
## the page is re-rendered by Ghostscript (pdfwrite) and, in an EndPage procedure, each lower-case tag is covered
## by a white box and the corresponding upper-case letter is drawn at the same position in Helvetica-Bold.
## Usage: uppercase_tags.pl in.pdf out.pdf minheight            (tags found with pdftotext -bbox)
##        uppercase_tags.pl in.pdf out.pdf glyph a:x0,y0,x1,y1 ... (tags given as glyph boxes in points, top-left origin)
use strict; use warnings;
my ($in, $out, $mode, @spec) = @ARGV;
my ($pw, $ph) = `pdfinfo "$in"` =~ /Page size:\s+([\d.]+) x ([\d.]+)/ or die "no page size\n";
my @tags;
if ($mode =~ /^[\d.]+$/) {
  for (split /\n/, `pdftotext -bbox "$in" -`) {
    next unless /xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([a-i])<\/word>/;
    my ($x0, $y0, $x1, $y1, $l) = ($1, $2, $3, $4, $5);
    next unless $y1 - $y0 >= $mode;
    my $fs = ($y1 - $y0) / 1.117;                    # font bounding box = ascent + descent (Arial/Helvetica)
    push @tags, [$l, $x0, $ph - $y1, $x1 - $x0, $y1 - $y0, $fs, $ph - $y1 + 0.212 * $fs];
  }
} else {
  for (@spec) { my ($l, $b) = split /:/; my ($x0, $y0, $x1, $y1) = split /,/, $b;
    my $fs = ($l =~ /[bdfhk]/ ? ($y1 - $y0) / 0.716 : ($y1 - $y0) / 0.519);   # ascender or x-height glyph
    push @tags, [$l, $x0, $ph - $y1, $x1 - $x0, $y1 - $y0, $fs, $ph - $y1]; }
}
die "no tags found\n" unless @tags;
my $ps = "<< /EndPage { exch pop 2 ne { gsave\n";
for my $t (@tags) {
  my ($l, $x, $yb, $w, $h, $fs, $base) = @$t;
  $ps .= sprintf("  1 setgray %.2f %.2f %.2f %.2f rectfill 0 setgray /Helvetica-Bold findfont %.2f scalefont setfont %.2f %.2f moveto (%s) show\n",
                 $x - 0.8, $yb - 0.8, $w + 1.6, $h + 1.6, $fs, $x, $base, uc $l);
  printf "%s -> %s at x=%.1f baseline=%.1f size=%.1f\n", $l, uc $l, $x, $base, $fs;
}
$ps .= "grestore true } { false } ifelse } bind >> setpagedevice\n";
open(my $o, '>', "$out.ps") or die; print $o $ps; close $o;
system("gs", "-q", "-o", $out, "-sDEVICE=pdfwrite", "-dDEVICEWIDTHPOINTS=$pw", "-dDEVICEHEIGHTPOINTS=$ph", "-dFIXEDMEDIA", "-dAutoRotatePages=/None", "$out.ps", "-f", $in) == 0 or die "gs failed\n";
print "wrote $out\n";
