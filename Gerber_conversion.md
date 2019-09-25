# How to Convert Cadence Allegro BRD to GDSII

1. Download full printed circuit board design from source
2. Unzip the design and navigate to the Gerber files in the electrical aspects
    * Gerber files may have any file extension but [the typical ones](https://www.pcbway.com/blog/help_center/Gerber_File_Extention_from_Different_Software.html) include .pcb, .art, .gbr or .ger, .top or .bot, and .gtl or .gbl
    * Excellon NC drill files likewise may have any file extension but the typical ones include .txt, .drl, .exc, .xln, .drd, and .tap
    * Unfortunately, most drill files use US customary units instead of SI units regardless of what measurement system the design software set, so scaling is typically needed as done below
3. Use [gerbv](http://gerbv.geda-project.org/), a free and open source tool for Gerber file viewing and limited editing that is part of the [gEDA project](http://www.geda-project.org/), to remove the mechanical drawing title block from each copper layer if it exists
4. Use [KLayout](https://github.com/KLayout/klayout) to import the Gerber files for each of the _N_ layers and the drill files for _plated through-holes_, all depths of _blind vias_, and all intervals of _buried vias_ into a single GDSII file (part of the import tool)
5. Scale the _N - 1_ via layers by 2.54 cm/in. and carefully reposition them so that they align perfectly with the _N_ metal layers
6. Save the resulting GDSII file to the working location with overlapping metal removed and discarding information that the older GDSII file format does not support but was imported from the Gerber files
7. Write a simulation input file in the working location with the stackup information from the design schematics
