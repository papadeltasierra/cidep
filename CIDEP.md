# Purpose
The port the CIDEP (Chemically Induced Dynamic Electron Polarisation) program from MS-DOS, FORTRAN-77 to a Jupyter Notebook and separately to Microsoft C#.

## AI Scope
-   AI may write, modify, delete files within the scope of the project but not outside it without permission
-   AI may access remote files to help with generation of new code
-   AI may move the files of, but must not modify, the existing FORTRAN-77 code
-   AI can restructure the layout of thye project directory so that it makes sense as a host for the three file sets
    -   Original FORTRAN-77
    -   Jupyter Notebook ported code
    -   C# ported code.
-   AI should be begin by creating a plan, indicating permissions that it requires, and asking for permission
    -   The intention is that once the user gives permission, the entire conversion process will happen without further user intervention
    -   Where AI makes decisions of its own, these should be noted for possible later modification
-   A new file, ./MIGRATION.md, should be created with sections indicating how each port was designed.

## Math Library
-   This program contains maths methods that can be replaced with standard open-source implementations.

## MS-DOS, FORTRAN-77 Overview
-   Source located in FORTRAN77/cidep directory
-   Some screen processing functions are in FORTRAN77/mode1
-	Single threaded
-	Simple textual menuing using DOS commands
-	Generates text input pages using ANSI escape codes
-	Generates graphical output in VGA or SVGA (800 x 600) mode using a library that interfaces to a pre-installed interupt driven library
-	Graphical cursor, with XOR colouring, is directed using keyboard keys and modifiers; there is no support for a mouse
-	A default set of configuration is provided for a specific radical pair
-	First phase generates a "stick diagram" of perfect spectra of signal strength to magnetic field strength
-	RPM (Radical Pair Mechanism) and TM (Triplet Mode) phases modifies the sinal strength as each mechanism would do
-   Spectrum shaping is added
-   Varying amounts of RPM to TM can be mixed to simulate an observed spectum (which is not accessed by this program)
-	Stores and reads data in a binary format

## New Requirements
-	Multi-threaded where possible and appropriate
-	Proper "Windows" menus
-	Input using standard "Windows" windows
-	Graphics using standard "Windows" graphical window
-	Cursor controlled by mouse but also accepting left-arrow and right-arrow when focussed to allow precise positioning where screen resolution is less than the resolution of the simulation
-	Store and read data in a JSON format not as binary
-   Ability to read, but not write, the original binary format

## C# vs Jupyter Notebook
-   Two implementations are to be created


