ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src data

dist_noinst_SCRIPTS = autogen.sh cleanup.sh

dist-hook:
	mkdir -p $(distdir)/resources/icons
	cp -rp $(srcdir)/resources/icons/*.png $(distdir)/resources/icons/
	cp -rp $(srcdir)/resources/icons/*.svg $(distdir)/resources/icons/
	mkdir -p $(distdir)/resources/icons/oxygen/16x16/actions
	cp -rp $(srcdir)/resources/icons/oxygen/AUTHORS $(distdir)/resources/icons/oxygen/
	cp -rp $(srcdir)/resources/icons/oxygen/COPYING $(distdir)/resources/icons/oxygen/
	cp -rp $(srcdir)/resources/icons/oxygen/index.theme $(distdir)/resources/icons/oxygen/
	cp -rp $(srcdir)/resources/icons/oxygen/16x16/actions/*.png $(distdir)/resources/icons/oxygen/16x16/actions/

include $(top_srcdir)/am/pkg-linux.am
include $(top_srcdir)/am/pkg-win.am