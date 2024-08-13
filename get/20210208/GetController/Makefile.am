include $(top_srcdir)/am/doxygen-ext.am

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src

EXTRA_DIST = doc/Doxyfile
dist_noinst_SCRIPTS = autogen.sh cleanup.sh

dist-hook:
	mkdir -p $(distdir)/res/icons
	cp -rp $(srcdir)/res/icons/*.png $(distdir)/res/icons/
	cp -rp $(srcdir)/res/icons/*.svg $(distdir)/res/icons/
	mkdir -p $(distdir)/res/icons/oxygen/16x16/actions
	cp -rp $(srcdir)/res/icons/oxygen/AUTHORS $(distdir)/res/icons/oxygen/
	cp -rp $(srcdir)/res/icons/oxygen/COPYING $(distdir)/res/icons/oxygen/
	cp -rp $(srcdir)/res/icons/oxygen/index.theme $(distdir)/res/icons/oxygen/
	cp -rp $(srcdir)/res/icons/oxygen/16x16/actions/*.png $(distdir)/res/icons/oxygen/16x16/actions/
