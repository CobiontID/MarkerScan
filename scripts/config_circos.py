from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-k", type=str, action='store', dest='karyo', metavar='karyo',help='define karyotype')
parser.add_argument("-d", type=str, action='store', dest='dat', metavar='dat',help='define data file')
parser.add_argument("-l", type=str, action='store', dest='link', metavar='link',help='define link file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

print('# circos.conf')

print('karyotype = '+results.karyo)

print('<ideogram>')
print('show_label     = yes')
print('label_with_tag = yes')
print('label_font     = light')
print('label_radius   = dims(ideogram,radius_outer) + 0.05r')
print('label_center   = yes')
print('label_size     = 48p')
print('label_color    = grey')
print('label_parallel = yes')
print('label_case     = upper')
#label_format   = eval(sprintf("chr%s",var(label)))

print('<spacing>')
print('default = 0.005r')
print('</spacing>')

print('radius    = 0.9r')
print('thickness = 20p')
print('fill      = yes')
print('</ideogram>')

print('<plots>')
print('        <plot>')
print('                show  = yes')
print('                type  = heatmap')
print('                file  = '+results.dat)
print('                r1    = 0.95r')
print('                r0    = 0.90r')
print('                color = red')
print('        </plot>')
print('</plots>')

print('<links>')
print('        <link>')
print('                 file = '+results.link)
print('               radius = 0.95r')
print('        bezier_radius = 0.1r')
print('                color = grey')
print('            thickness = 3')
#fill_color = red
print('        </link>')
print('</links>')
################################################################
################################################################
# The remaining content is standard and required. It is imported
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files,
# look in etc/ in the Circos distribution.

print('<image>')
# Included from Circos distribution.
print('<<include etc/image.conf>>')
print('</image>')

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
print('<<include etc/colors_fonts_patterns.conf>>')

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
print('<<include etc/housekeeping.conf>>')