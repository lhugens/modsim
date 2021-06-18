MGL_LIBS=-L/usr/X11/lib -L/usr/local/lib -L/usr/X11R6/lib -L/usr/lib -lm -lXmu -lXi  -lX11 -L/sw/lib -lGLU -lGL -lglut -lX11 -lXext -L/sw/lib  -lpng
EXTRA_OPTIONS_FOR_MAC=-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib: /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib

MGL_INCLUDE= -I$(MD)/include -I/opt/local/include -I/usr/X11R6/include -I/usr/local/include -UMGL_USELIST -I/sw/include

molgl: molgl.o #mglgui.o
	gcc -lpng -O3 -o molgl molgl.o $(MGL_LIBS) #$(EXTRA_OPTIONS_FOR_MAC)

molgl.o: molgl.c molgl.h
	gcc -O3 -c -Wall -UMGL_MACOS $(MGL_INCLUDE) -o molgl.o molgl.c
#USARE QUESTO PER MAC OSX: gcc -O0 -g -c -Wall -DMGL_MACOS $(MGL_INCLUDE) -o molgl.o molgl.c

#mglgui.o : mglgui.C mglgui.H
#	g++ -c -Wall $(MGL_INCLUDE) -o mglgui.o mglgui.C

clean:
	-rm *.o *~

realclean:
	-rm molgl
	-rm *.o *~
