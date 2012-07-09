ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = erne-bs5 erne-create erne-filter erne-map
noinst_LIBRARIES = liberne.a liberne-sam.a
 
if DEBUG_ENABLED
bin_PROGRAMS += erne-bs5-gdb erne-create-gdb erne-filter-gdb erne-map-gdb
noinst_LIBRARIES += liberne-gdb.a liberne-sam-gdb.a
endif

if HAVE_MPI
MPI_CFLAGS = `$(MPICC) --showme:compile`
MPI_CXXFLAGS = `$(MPICXX) --showme:compile`
MPI_LFLAGS = `$(MPICXX) --showme:link`

bin_PROGRAMS += erne-dmap erne-dcreate
noinst_LIBRARIES += libderne.a libderne-sam.a
if DEBUG_ENABLED
bin_PROGRAMS += erne-dmap-gdb erne-dcreate-gdb
noinst_LIBRARIES += libderne-gdb.a libderne-sam-gdb.a
endif
endif

if CYGWIN_TARGET
install-exec-hook:
	cp -p /usr/bin/cygwin1.dll $(DESTDIR)$(bindir) 
endif

ld_generic_flags = -L$(BOOST_LDPATH) -R$(BOOST_LDPATH) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_REGEX_LIBS) -lz -lbz2

##### LIBS ######

liberne_a_SOURCES = \
	options/Options.cpp \
	data_structures/GlobalToLocal.cpp \
	data_structures/Hash.cpp \
	data_structures/Mask.cpp \
	data_structures/SetZ.cpp \
	data_structures/MismatchCombination.cpp \
	errors/Data_Exception.cpp \
	errors/Data_Not_Found.cpp \
	errors/File_Not_Found.cpp \
	errors/Generic_Exception.cpp \
	errors/Incorrect_Format.cpp \
	io/Auto_Unzip.cpp \
	io/Fasta.cpp \
	io/Samfile.cpp \
	io/read_functions.cpp
liberne_a_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS)
liberne_gdb_a_SOURCES = $(liberne_a_SOURCES)
liberne_gdb_a_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)

libderne_a_SOURCES = $(liberne_a_SOURCES)
libderne_a_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
libderne_gdb_a_SOURCES = $(libderne_a_SOURCES)
libderne_gdb_a_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)

liberne_sam_a_SOURCES = \
	samtools/bam_aux.c \
	samtools/bam.c \
	samtools/bam_import.c \
	samtools/bam_index.c \
	samtools/bam_pileup.c \
	samtools/bgzf.c \
	samtools/faidx.c \
	samtools/knetfile.c \
	samtools/kstring.c \
	samtools/razf.c \
	samtools/sam.c \
	samtools/sam_header.c \
	samobjects/BAM1.cpp \
	samobjects/SAM_IO_Error.cpp \
	samobjects/SAM_istream.cpp \
	samobjects/SAM_ostream.cpp \
	samobjects/SAM_search.cpp \
	samobjects/SAM.cpp \
	samobjects/Temporary_File.cpp
liberne_sam_a_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS)
liberne_sam_CFLAGS = -O3 -fno-stack-protector
liberne_sam_gdb_a_SOURCES = $(liberne_sam_a_SOURCES)
liberne_sam_gdb_a_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)
liberne_sam_gdb_CFLAGS =  -O0 -g3 -fno-stack-protector

libderne_sam_a_SOURCES = $(liberne_sam_a_SOURCES)
libderne_sam_a_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
libderne_sam_CFLAGS = -O3 -fno-stack-protector $(MPI_CFLAGS)
libderne_sam_gdb_a_SOURCES = $(libderne_sam_a_SOURCES)
libderne_sam_gdb_a_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
libderne_sam_gdb_CFLAGS =  -O0 -g3 -fno-stack-protector $(MPI_CFLAGS)

##### erne-bs5 #####
erne_bs5_SOURCES = erne-bs5.cpp \
	options/Options_BS5.cpp \
	modules/Module_MAP.cpp \
	modules/Module_BS5.cpp
erne_bs5_DEPENDENCIES = liberne.a liberne-sam.a
erne_bs5_CXXFLAGS = -O3  -fno-stack-protector $(BOOST_CPPFLAGS)
erne_bs5_LDFLAGS = $(erne_bs5_DEPENDENCIES) $(ld_generic_flags)
if STATIC_ENABLED
erne_bs5_LDFLAGS += -all-static
endif
erne_bs5_LIBS = $(erne_bs5_LDFLAGS)

erne_bs5_gdb_DEPENDENCIES = liberne-gdb.a liberne-sam-gdb.a
erne_bs5_gdb_SOURCES = $(erne_bs5_SOURCES)
erne_bs5_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)
erne_bs5_gdb_LDFLAGS = $(erne_bs5_gdb_DEPENDENCIES) $(erne_bs5_LDFLAGS)
erne_bs5_gdb_LIBS = $(erne_bs5_LDFLAGS)

##### erne-create #####
erne_create_SOURCES = erne-create.cpp \
	options/Options_CREATE.cpp \
	modules/Module_CREATE.cpp
	
erne_create_DEPENDENCIES = liberne.a
erne_create_CXXFLAGS = -O3  -fno-stack-protector $(BOOST_CPPFLAGS)
erne_create_LDFLAGS = $(erne_create_DEPENDENCIES) $(ld_generic_flags)
if STATIC_ENABLED
erne_create_LDFLAGS += -all-static
endif
erne_create_LIBS = $(erne_create_LDFLAGS)

erne_create_gdb_DEPENDENCIES = liberne-gdb.a
erne_create_gdb_SOURCES = $(erne_create_SOURCES)
erne_create_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)
erne_create_gdb_LDFLAGS = $(erne_create_gdb_DEPENDENCIES) $(erne_create_LDFLAGS)
erne_create_gdb_LIBS = $(erne_create_LDFLAGS)

##### erne-filter #####
erne_filter_SOURCES = erne-filter.cpp \
	options/Options_FILTER.cpp \
	modules/Module_FILTER.cpp

erne_filter_DEPENDENCIES = liberne.a
erne_filter_CXXFLAGS = -O3  -fno-stack-protector $(BOOST_CPPFLAGS)
erne_filter_LDFLAGS = $(erne_filter_DEPENDENCIES) $(ld_generic_flags)
if STATIC_ENABLED
erne_filter_LDFLAGS += -all-static
endif
erne_filter_LIBS = $(erne_filter_LDFLAGS)

erne_filter_gdb_DEPENDENCIES = liberne-gdb.a
erne_filter_gdb_SOURCES = $(erne_filter_SOURCES)
erne_filter_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)
erne_filter_gdb_LDFLAGS = $(erne_filter_gdb_DEPENDENCIES) $(erne_filter_LDFLAGS)
erne_filter_gdb_LIBS = $(erne_filter_LDFLAGS)

##### erne-map #####
erne_map_SOURCES = erne-map.cpp \
	options/Options_MAP.cpp \
	modules/Module_MAP.cpp

erne_map_DEPENDENCIES = liberne.a liberne-sam.a
erne_map_CXXFLAGS = -O3  -fno-stack-protector $(BOOST_CPPFLAGS)
erne_map_LDFLAGS = $(erne_map_DEPENDENCIES) $(ld_generic_flags)
if STATIC_ENABLED
erne_map_LDFLAGS += -all-static
endif
erne_map_LIBS = $(erne_map_LDFLAGS)

erne_map_gdb_DEPENDENCIES = liberne-gdb.a liberne-sam-gdb.a
erne_map_gdb_SOURCES = $(erne_map_SOURCES)
erne_map_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS)
erne_map_gdb_LDFLAGS = $(erne_map_gdb_DEPENDENCIES) $(erne_map_LDFLAGS)
erne_map_gdb_LIBS = $(erne_map_LDFLAGS)

if HAVE_MPI

##### erne-dcreate #####
erne_dcreate_SOURCES = erne-dcreate.cpp \
	options/Options_DCREATE.cpp \
	modules/Module_DCREATE.cpp

erne_dcreate_DEPENDENCIES = libderne.a
erne_dcreate_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
erne_dcreate_LDFLAGS = $(erne_dcreate_DEPENDENCIES) $(ld_generic_flags) $(MPILIBS) $(MPI_LFLAGS)
erne_dcreate_LIBS = $(erne_dcreate_LDFLAGS)

erne_dcreate_gdb_DEPENDENCIES = libderne-gdb.a
erne_dcreate_gdb_SOURCES = $(erne_dcreate_SOURCES)
erne_dcreate_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
erne_dcreate_gdb_LDFLAGS = $(erne_dcreate_gdb_DEPENDENCIES) $(ld_generic_flags) $(MPILIBS) $(MPI_LFLAGS)
erne_dcreate_gdb_LIBS = $(erne_dcreate_gdb_LDFLAGS)

##### erne-dmap #####
erne_dmap_SOURCES = erne-dmap.cpp \
	io/Hash_To_Samfile.cpp \
	options/Options_DMAP.cpp \
	modules/Module_MAP.cpp \
	modules/Module_DMAP.cpp

erne_dmap_DEPENDENCIES = libderne.a libderne-sam.a
erne_dmap_CXXFLAGS = -O3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
erne_dmap_LDFLAGS = $(erne_dmap_DEPENDENCIES) $(ld_generic_flags) $(MPILIBS) $(MPI_LFLAGS)
erne_dmap_LIBS = $(erne_dmap_LDFLAGS)

erne_dmap_gdb_DEPENDENCIES = libderne-gdb.a libderne-sam-gdb.a
erne_dmap_gdb_SOURCES = $(erne_dmap_SOURCES)
erne_dmap_gdb_CXXFLAGS = -O0 -g3 -fno-stack-protector $(BOOST_CPPFLAGS) -DHAVE_MPI $(MPI_CXXFLAGS)
erne_dmap_gdb_LDFLAGS = $(erne_dmap_gdb_DEPENDENCIES) $(ld_generic_flags) $(MPILIBS) $(MPI_LFLAGS)
erne_dmap_gdb_LIBS = $(erne_dmap_gdb_LDFLAGS)

endif
