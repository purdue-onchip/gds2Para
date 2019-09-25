# How to Convert Cadence Allegro BRD to GDSII

1. Locate the Cadence Allegro PCB file (.brd extension)
2. Run Cadence DB Doctor (separately installed program in the Windows Start Menu) to update the file to the installed version of Cadence Allegro to be compatible with Cadence PCB Editor
3. Select the input design using the ellipsis to open a dialog box with a file manager
4. Select the output design name using the ellipsis to open a dialog box with a file manager (a different name is recommended to preserve the original file)
5. Check "Update all DRC (including Batch)", "Check shape outlines", and "Regenerate Xnets" before pressing the "Check" button
6. Close out of Cadence DB Doctor if the logs say the file conversion is successful
7. Run Cadence PCB Editor to open a dialog box for Cadence Allegro
8. Select the product "Allegro PCB Designer" with none of the options checked, and then press "OK"
9. Right click on the toolbars at the top of Allegro PCB Designer and ensure that "Visibility", "Command", and "View" above the divider as well as "File"', "Edit", "View" (again), and "Manufacture" below the divider are enabled so that these toolbars appear
10. Use the menu option "File" &rightarrow; "Open" to load the output design name used with Cadence DB Doctor
11. Click "On" for global visibility on the "Visibility" toolbar on the right side of the screen because all features must be visible in order to export to a different file format
12. Use the menu option "Setup" &rightarrow; "Cross-section..." to open the Cross Section Editor dialog box and record the layer properties (the "Export" menu lets you export the information to an XML file if desired)
13. Close the Cross Section Editor dialog box when finished
14. Use the menu option "Manufacture" &rightarrow; **"Stream Out"** (note the name) to open the Stream Out dialog box which is the GDSII file export wizard
15. Select the output GDSII file name using the ellipsis to open a dialog box with a file manager
16. Select a layer conversion file name using the ellipsis to open a dialog box with a file manager (this file name with the .cnv extension should not already be in used at the target location)
17. Click "Edit..." below the layer conversion file name to open the Stream Out Edit Layer Conversion File dialog box
18. Use the class filter to bring up the following classes and subclasses for conversion, checking the "Selected" box in the table for each subclass listed before adjusting the map selected items "Layer" number in its box starting from 1 and clicking "Map"
    * _Etch_
        1. Every subclass corresponding to the layers given in the Cross Section Editor dialog box (recommended to use only odd numbers starting from 1 in order if feasible)
        2. Skip subclasses not corresponding to layers (e.g., _WIRE_)
    * _Via\_Class_
        1. Every subclass corresponding to the layers given in the Cross Section Editor dialog box (recommended to use only even numbers starting from 2 in order if feasible)
        2. Skip subclasses not corresponding to layers (e.g., _FILMMASKBOTTOM_, _PASTEMASK\_BOTTOM_, and _SOLDERMASK\_TOP_)
    * _Board\_Geometry_
        1. Only the _OUTLINE_ subclass (useful for the future of the simulation input file)
19. Back in the Stream Out dialog box, select the endcap style for GDSII path elements in the drop-down menu (e.g., _Flush_)
20. Check "Add pin text to top level in hierarchy" and leave the other default options in place
21. Click the "Export" button and close everything without saving in Allegro PCB Editor if the logs say the export was successful
22. Save the resulting GDSII file to the working location, manipulating the file with [KLayout layer Boolean operations](https://www.klayout.de/doc/manual/layer_boolean.html) if there are issues in appearance
23. Write a simulation input file in the working location with the stackup information from the Cross Section Editor dialog box and mapping of the stream out layer conversion file (ASCII file with .cnv extension)
