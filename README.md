# public_PathReporter
Creating Value from Free Text and Semi-structure Pathology Reports (AGPL v3.0)
=============================================================================================
    AGPL v3.0 NOTICE
    PLEASE INCLUDE AT START OF ALL SOURCE CODE FILES PERTAINING TO THIS PROGRAM

    PROGRAM TITLE: PathReporter
    Copyright (C) 2022  Hansen Lam, M.D., Freddy Nguyen, M.D., Ph.D., Karl Eichorn, Dan Wild
    Disclaimer: This software is intended for research purposes only. This software is not
    for clinical use. The authors and contributors shall not be liable for any damages or
    losses resulting from the use of this software. The authors and contributors shall not 
    be liable for any actions or decisions resulting rom the use of this software.
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


OBJECTIVE:
This program processes semi-structured pathology reports in XML format to create tabulated data and
facilitate data collection from free-text reports. 

METHODS:
This program implements an object oriented approach, using Python v3.x. Dependencies include but are
not limited to:
pandas, numpy, openpyxl, re, string, math, lxml, ujson, collections
