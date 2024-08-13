
getbenchdir = $(datadir)/get-bench
getbench_DATA = 


hardware_files = config/hardwareDescription_fullCoBoStandAlone.xcfg \
 config/hardwareDescription_Mutant.xcfg \
 config/describe-mono-node.xcfg config/log4j.properties \
 config/configure-pedestals.xcfg \
 config/configure-ExternalTrigger-2CoBo_1Mutant.xcfg
cfgdir = $(getbenchdir)/config
cfg_DATA = $(hardware_files)

format_files = \
	format/CoboFormats.xcfg \
	format/CoboFormats-Rev-0.xcfg \
	format/CoboFormats-Rev-1.xcfg \
	format/CoboFormats-Rev-2.xcfg \
	format/CoboFormats-Rev-3.xcfg \
	format/CoboFormats-Rev-4.xcfg \
	format/CoboFormats-Rev-5.xcfg \
	format/CoboFormats-Rev-5-Compact.xcfg \
	format/MergedDataFormats-ByEventTime-Rev-1.xcfg \
	format/MergedDataFormats-ByEventId-Rev-1.xcfg
formatdir = $(getbenchdir)/format
format_DATA = $(format_files)

ecc_files = \
	scripts/ice/adc.ecc \
	scripts/ice/adcCal.ecc \
	scripts/ice/adcTest.ecc \
	scripts/ice/aget.ecc \
	scripts/ice/agetBaseLine.ecc \
	scripts/ice/agetBaseLineFIFOBug.ecc \
	scripts/ice/agetCal.ecc \
	scripts/ice/asad.ecc \
	scripts/ice/asadCal.ecc \
	scripts/ice/asadInjectCharge.ecc \
	scripts/ice/beforeStart.ecc \
	scripts/ice/data.ecc \
	scripts/ice/initEverything.ecc \
	scripts/ice/lcm_phase0.ecc \
	scripts/ice/machine.ecc \
	scripts/ice/pipe.ecc \
	scripts/ice/powerUp.ecc \
	scripts/ice/sequencer.ecc \
	scripts/ice/switchToFastControl.ecc \
	scripts/ice/switchToSlowControl.ecc \
	scripts/ice/cobo_pll_100MHz_CoBo0.2_StandaloneConfig.ecc \
	scripts/ice/cobo_pll_100MHz_CoBo1.0_StandaloneConfig.ecc \
	scripts/ice/cobo_pll_100MHz_MutantConfig.ecc \
	scripts/ice/asad-adc-allones.ecc \
	scripts/ice/asad-adc-allzeros.ecc \
	scripts/ice/asad-adc-custom1000.ecc \
	scripts/ice/asad-adc-deskew.ecc \
	scripts/ice/asad-adc-normal.ecc
	
scripts_eccdir = $(getbenchdir)/scripts/ice
scripts_ecc_DATA = $(ecc_files)

soap_scripts = scripts/soap/soap_client_run_control_scripts.sh
scripts_soapdir = $(getbenchdir)/scripts/soap
scripts_soap_DATA = $(soap_scripts)


EXTRA_DIST = $(hardware_files) $(format_files) $(ecc_files) $(soap_scripts)
