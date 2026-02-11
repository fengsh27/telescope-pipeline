#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  24 10:29:10 2022

@author: Rosario Distefano
"""

""" Resources """
from lib.cmd import CmdParser
from lib.worker_patch import PatchedWorker

if __name__ == '__main__':
            
    cmd_obj: CmdParser = CmdParser()
    wk_obj: PatchedWorker = PatchedWorker(cmd_obj)
    wk_obj.run()
