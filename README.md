potracej
========

Java port of potrace (main parts).  Original code: http://potrace.sourceforge.net
This particular branch to port was extracted from current inkscape.

Contains decomposer and curve optimizer, no backends.

Ported from inkscape, java-style path translator.

Also simple UI to test it is working.

Also, ant makefile.

WARNING! NON-DEFAULT path is not tested.

---------------------------

Contains one less then optimal part in decomposer, poorly converted to java due to lack of pointers. 
Was too tired to sort it out, decomposer is 25% of all CPU work, maybe could be reduced to 10; enough;


