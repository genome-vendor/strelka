
STRELKA_DIR := $(CURDIR)/strelka
export REDIST_DIR := $(CURDIR)/redist
export BOOST_ID := boost_1_44_0_subset
export CODEMIN_ID := CodeMin-1.0.2
export SAMTOOLS_ID := samtools-0.1.14_no_tview
export TABIX_ID := tabix-0.2.5_new_libz_order
export VCFTOOLS_ID := vcftools_0.1.6

LAUNCH_SCRIPT := configureStrelkaWorkflow.pl

all:
	$(MAKE) -C $(REDIST_DIR) && \
	$(MAKE) -C $(STRELKA_DIR) && \
        cp $(STRELKA_DIR)/src/perl/bin/$(LAUNCH_SCRIPT) $(LAUNCH_SCRIPT)

clean: srcclean
	$(MAKE) -C $(REDIST_DIR) clean

# developer targets:
#

# clean only strelka directories but leave redist built:
srcclean:
	$(MAKE) -C $(STRELKA_DIR) clean
	rm -f $(LAUNCH_SCRIPT)

etags:
	cd $(STRELKA_DIR) && find . -type f -iname "*.[ch]" -or -iname "*.cpp" -or -iname "*.hh" | sed "s/\.\///" | etags -	
