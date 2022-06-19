from pymol import cmd
import argparse
import os


fps = 20
mode = 0

if mode == 0:
    for i in range(1,101):
        cmd.frame(f'{i}')
        cmd.ray(1000,1000)
        cmd.png(f'{i:05d}.png')

delay = 100//fps
os.system(f'convert -delay {delay} -dispose Background *.png anim.gif')