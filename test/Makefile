
SUBDIRS :=$(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN :=$(addsuffix .clean,$(SUBDIRS))
.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS)
	make/check

# this is cleaning all subdirectories
clean: $(SUBDIRSCLEAN)

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean
$(SUBDIRS):
	$(MAKE) -C $@






