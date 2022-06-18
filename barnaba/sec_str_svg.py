#!/usr/bin/env python

import numpy as np
from matplotlib.pylab import *
import barnaba.sec_str_constants as secon

n_colors = 128
cmap = matplotlib.cm.get_cmap('gnuplot_r', n_colors) 
def hexcolor(value):
    rgb = cmap(int(value*n_colors)) 
    return matplotlib.colors.rgb2hex(rgb)

def draw_structure(threshold, pos, pairs, ann_list, chi_conf, sequence, dimensions, output_ids):
    binary = False
    syn = np.fromstring(np.binary_repr(chi_conf), dtype='S1').astype(int)[::-1]
    n = len(pos)
    output_svg = ""

    output_svg += draw_header(dimensions)
    output_foreground = ""
    output_background = ""
    for i in range(n-1):
        search = (i, i+1)
        if len(np.where((pairs == search).all(axis=1))[0]) == 0:
            xy1 = pos[i]
            xy2 = pos[i+1]
            output_svg += draw_dummy(xy1, xy2)
    for i, pair in enumerate(pairs):
        n = sum([1 for key, value in ann_list.items() if pair[0] in key and pair[1] in key and value >= threshold])
        if n == 0:
            continue
        if n > 1:
            k = 0
            d = (n-1)*1.6
        elif n > 2:
            print("MORE than 2 annotations for pair %s!!" % pair)
        for key, value in ann_list.items():
            if not pair[0] in key or not pair[1] in key:
                continue
            ann = key[2]
            if value >= threshold:
                color = hexcolor(value)
                r1 = pair[0] 
                r2 = pair[1]
                x1 = pos[r1][0]
                y1 = pos[r1][1]
                x2 = pos[r2][0]
                y2 = pos[r2][1]
                if n > 1 :
                    rot = np.arctan2((x1-x2), (y2-y1)) + np.pi*.5
                    dx = d * np.sin(rot)
                    dy = d * np.cos(rot)
                    if k == 0:
                        x1 += dx
                        y1 -= dy
                        x2 += dx
                        y2 -= dy
                    elif (n == 2 and k == 1) or (n == 3 and k == 2):
                        x1 -= dx
                        y1 += dy
                        x2 -= dx
                        y2 += dy
                    k += 1
                if ann in secon.list_stackings:
                    output_foreground += draw_stack([x1, y1], [x2, y2], ann, color)
                else:
                    if ann == "GUc":
                        t_ann = "WWc"
                    elif ann == "WCc" and ((sequence[r1][0] == "C" and sequence[r2][0] == "G") or
                            (sequence[r1][0] == "G" and sequence[r2][0] == "C")):
                        t_ann = "=c"
                    elif ann == "WCc" and ((sequence[r1][0] == "A" and sequence[r2][0] == "U") or
                            (sequence[r1][0] == "U" and sequence[r2][0] == "A")):
                        t_ann = "-c"
                    else:
                        t_ann = ann
                    if t_ann == "WCc":
                        print("Found WC annotation for pair ", r1, r2, sequence[r1], sequence[r2])
                    if binary and ann in secon.list_wc_pairs:
                        color = "#ff0000"
                    if t_ann in ["-c", "=c"] or t_ann[2] == "t":
                        output_foreground += draw_basepair([x1, y1], [x2, y2], t_ann, color)
                    else:
                        output_background += draw_basepair([x1, y1], [x2, y2], t_ann, color)

    output_svg += output_background
    output_svg += output_foreground
             
    for i, p in enumerate(pos):
        try: syn[i]
        except: background = "#ffffff"
        else:
            if syn[i] == 1:
                background = "#00ff00"
            else:    
                background = "#ffffff"
        if output_ids:
           # output_svg += draw_res(p, sequence[i][1], background)
            output_svg += draw_res(p, i, background)
        else:    
            output_svg += draw_res(p, sequence[i][0], background)    

    output_svg += draw_footer()
    return output_svg


def draw_header(length):
    output_header = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<svg version="1.1"
   width="%.0fpx"
   height="%.0fpx"
xmlns="http://www.w3.org/2000/svg">
    <rect width="100%%" height="100%%" fill="white"/>
''' % (length, length)
    return output_header


def draw_footer():
    return "\n</svg>"

def draw_number(pos, res_shift=1):
    output_number = ""
    output_number += '''
<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="rgb(25%%, 25%%, 25%%)" stroke-width="1.0" />''' % (pos[0][0], pos[0][1]+12.5, pos[0][0], pos[0][1]+7.5)
    output_number += '''
<text x="%.4f" y="%.4f" text-anchor="middle" font-family="Helvetica" font-size="7.5" fill="rgb(25%%, 25%%, 25%%)" >%d</text>''' % (pos[0][0], pos[0][1]+23, res_shift)
    output_number += '''
<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="rgb(25%%, 25%%, 25%%)" stroke-width="1.0" />''' % (pos[-1][0], pos[-1][1]+12.5, pos[-1][0], pos[-1][1]+7.5)
    output_number += '''
<text x="%.4f" y="%.4f" text-anchor="middle" font-family="Helvetica" font-size="7.5" fill="rgb(25%%, 25%%, 25%%)" >%d</text>''' % (pos[-1][0], pos[-1][1]+23, len(pos)-1+res_shift)
    output_number += '''
<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="rgb(25%%, 25%%, 25%%)" stroke-width="1.0" />''' % (pos[9][0]-12.5, pos[9][1], pos[9][0]-7.5, pos[9][1])
    output_number += '''
<text x="%.4f" y="%.4f" text-anchor="middle" font-family="Helvetica" font-size="7.5" fill="rgb(25%%, 25%%, 25%%)" >%d</text>''' % (pos[9][0]-20, pos[9][1]+3, 9+res_shift)
    output_number += '''
<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="rgb(25%%, 25%%, 25%%)" stroke-width="1.0" />''' % (pos[19][0]+12.5, pos[19][1], pos[19][0]+7.5, pos[19][1])
    output_number += '''
<text x="%.4f" y="%.4f" text-anchor="middle" font-family="Helvetica" font-size="7.5" fill="rgb(25%%, 25%%, 25%%)" >%d</text>''' % (pos[19][0]+20, pos[19][1]+3, 19+res_shift)
    return output_number

def draw_res(xy, letter, background="#ffffff"):
    x = xy[0]
    y = xy[1]
    output_res = '''
<circle cx="%.4f" cy="%.4f" r="5.0" stroke="#595959" stroke-width="1.0" fill="%s"/>
<text x="%.4f" y="%.4f" text-anchor="middle" font-family="Helvetica" font-size="7.5" fill="#000000" >%s</text>''' % (x, y, background, x, y+2.7, letter)
    return output_res

def draw_dummy(xy1, xy2):
    x1 = xy1[0]
    y1 = xy1[1]
    x2 = xy2[0]
    y2 = xy2[1]
    color = "#595959"
    output = '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s"
     style="stroke-width:2;stroke-miterlimit:4;stroke-dasharray:1,1;stroke-dashoffset:0" />''' % (x1, y1, x2, y2, color)
    return output 

def draw_stack(xy1, xy2, _type, color):
    x1 = xy1[0]
    y1 = xy1[1]
    x2 = xy2[0]
    y2 = xy2[1]
    x_c = .5 * (x1 + x2)
    y_c = .5 * (y1 + y2)
    rot = np.arctan2((x1-x2), (y2-y1))/np.pi*180+90
    output = ""
    output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
    if _type == "<>":
        output += '''
  <rect
     style="fill:#ffffff;stroke:none"
     width="5.5"
     height="2"
     x="%.4f"
     y="%.4f" 
     transform="rotate(%.4f %.4f %.4f)"    /> ''' % (x_c-2.75, y_c-1, rot, x_c, y_c)
        output += '''
  <path
     style="fill:none;stroke:%s;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
    d="m %.4f,%.4f 2,2 -2,2 m -1.5,-4 -2,2 2,2" transform="rotate(%.4f %.4f %.4f)"/>''' % (color, x_c+.75, y_c-2, rot, x_c, y_c)
    if _type == ">>" or _type == "<<":
        if _type == "<<": rot += 180
        output += '''
  <rect
     style="fill:#ffffff;fill-opacity:1;fill-rule:nonzero;stroke:none"
     width="2"
     height="2"
     x="%.4f"
     y="%.4f"
     transform="rotate(%.4f %.4f %.4f)" /> ''' % (x_c-.5, y_c-1, rot, x_c, y_c) 
        output += '''
  <path
     style="fill:none;stroke:%s;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
     d="m %.4f,%.4f 2,2 -2,2 m -2.370121,-4 2,2 -2,2"
     transform="rotate(%.4f %.4f %.4f)"/>''' % (color, x_c, y_c-2, rot, x_c, y_c)
    if _type == "><":
        output += '''
  <rect
     style="fill:#ffffff;fill-opacity:1;fill-rule:nonzero;stroke:none"
     width="2"
     height="2"
     x="%.4f"
     y="%.4f" 
     transform="rotate(%.4f %.4f %.4f)" /> ''' % (x_c-1, y_c-1, rot, x_c, y_c) 
        output += '''
  <path
     style="fill:none;stroke:%s;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
     d="m %.4f,%.4f -2,2 2,2 m -6,-4 2,2 -2,2"
     transform="rotate(%.4f %.4f %.4f)"/>''' % (color, x_c+3, y_c-2, rot, x_c, y_c)
    return output


def draw_basepair(xy1, xy2, _fulltype, color):
    x1 = xy1[0]
    y1 = xy1[1]
    x2 = xy2[0]
    y2 = xy2[1]
    _type = _fulltype[:-1]
    ct_type = _fulltype[-1]
    if ct_type == "t":
    #    print "TRANS BASEPAIR! ", _fulltype
        fill = "#ffffff"
    else:
        fill = color
    x_c = .5 * (x1 + x2)
    y_c = .5 * (y1 + y2)
    rot = np.arctan2((x1-x2), (y2-y1))/np.pi*180+90
    rot_rad = rot/180*np.pi
#    print rot
    output = ""
    if _type == "WW":
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <circle cx="%.4f" cy="%.4f" r="2.25" stroke="%s" stroke-width="1.0" fill="%s"/>''' % (x_c, y_c, color, fill)
    elif _type == "-":
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
    elif _type == "=":
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1-1.25*np.sin(rot_rad), y1+1.25*np.cos(rot_rad), x2-1.25*np.sin(rot_rad), y2+1.25*np.cos(rot_rad), color)
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1+1.25*np.sin(rot_rad), y1-1.25*np.cos(rot_rad), x2+1.25*np.sin(rot_rad), y2-1.25*np.cos(rot_rad), color)
    elif _type == "SW" or _type == "WS":
        if _type == "WS": rot += 180
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <circle cx="%.4f" cy="%.4f" r="2.25" stroke="%s" stroke-width="1.0" fill="%s"
    transform="rotate(%.4f %.4f %.4f)" />''' % (x_c+3.5, y_c, color, fill, rot, x_c, y_c)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s;stroke-width:1.0" fill="%s" 
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c-1.75, y_c+2.25, x_c-1.75, y_c-2.25, x_c-5.75, y_c, color, fill, rot, x_c, y_c)

    elif _type == "HW" or _type == "WH":
        if _type == "WH": rot += 180
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <circle cx="%.4f" cy="%.4f" r="2.25" stroke="%s" stroke-width="1.0" fill="%s"
    transform="rotate(%.4f %.4f %.4f)" />''' % (x_c+3.75, y_c, color, fill, rot, x_c, y_c)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s; stroke-width:1.0" fill="%s" 
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c-6.0, y_c+2.25, x_c-6.0, y_c-2.25, x_c-1.5, y_c-2.25, x_c-1.5, y_c+2.25, color, fill, rot, x_c, y_c)
    elif _type == "HS" or _type == "SH":
        if _type == "SH": rot += 180
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s;stroke-width:1.0" fill="%s" 
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c-5.75, y_c+2.25, x_c-5.75, y_c-2.25, x_c-1.25, y_c-2.25, x_c-1.25, y_c+2.25, color, fill, rot, x_c, y_c)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s;stroke-width:1.0" fill="%s"
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c+1.75, y_c+2.25, x_c+1.75, y_c-2.25, x_c+5.75, y_c, color, fill, rot, x_c, y_c)
    elif _type == "HH":
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s;stroke-width:1.0" fill="%s" 
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c-2.25, y_c+2.25, x_c-2.25, y_c-2.25, x_c+2.25, y_c-2.25, x_c+2.25, y_c+2.25, color, fill, rot, x_c, y_c)
    elif _type == "SS":
        output += '''
  <line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="1.0" />''' % (x1, y1, x2, y2, color)
        output += '''
  <path d="M %.4f %.4f L %.4f %.4f L %.4f %.4f z" style="stroke:%s;stroke-width:1.0" fill="%s"
    transform="rotate(%.4f %.4f %.4f)"/>''' % (x_c-2, y_c+2.25, x_c-2, y_c-2.25, x_c+2, y_c, color, fill, rot, x_c, y_c)
    
    else:
        print("type %s unknown" % _type)
    return output
