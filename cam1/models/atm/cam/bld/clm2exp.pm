#
#	clm2exp.pm			Erik Kluzek
#
#	Perl module to create a namelist for CLM2.
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package clm2exp;

use atmlndnl;
@clm2exp::ISA = qw(atmlndnl  namelist);

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $CAM_config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'out'};
  my $printlev = $$optsref{'printlev'};

  my $default_vals = {};

  my $self = $class->SUPER::new( "clmexp", \%main::CLMEXP, $interactive, $file, 
                                 "DefaultCLMEXPNamelist.xml", $default_vals,
                                 $CAM_config, $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  $self->{DYNAMICS}   = $CAM_config->cfg("DYNAMICS");       # Dynamics to use
  $self->{RESOLUTION} = $CAM_config->cfg("RESOLUTION");     # horizontal resolution
  $self->{OCEANMODEL} = $CAM_config->cfg("OCEANMODEL");     # Ocean using with CAM (dom or som)

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the CLM2 namelist variables.

  my ($self, %settings) = @_;

  my $runtype = $settings{RUNTYPE};
  my $class = ref($self);
  my $nm = "$class\:\:set_default_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $default_vals = $self->{'default_vals'};
  my $opt;
  my $raw_files_spec;

  # Get the default values from the XML file
  $self->get_default_values( %settings );

  # Check that "nrevsn" is set if this is a branch simulation
  if ($runtype eq 'branch' and !defined($NLref->{'nrevsn'})) {
      if ( $self->do_interactive ) {
	  print "Enter absolute pathname for CLM2 master restart file from which to branch: ";
	  $opt = <>; chomp $opt;
	  $NLref->{'nrevsn'} = namelist::quote_string($opt);
      } else {
	  die "ERROR: The CLM2 master restart file must be specified for a branch\n".
	  "       run.  Set the namelist variable NREVSN to the absolute\n".
	  "       pathname for this dataset.\n".
	  "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option.\n";
      }
  }

  # Root directory for default initial and boundary datasets
  my $rootdir;
  if (defined($optsref->{'csmdata'})) {
      $rootdir = $optsref->{'csmdata'};
  } elsif (defined $ENV{'CSMDATA'}) {
      $rootdir = $ENV{'CSMDATA'};
  } else {
      $rootdir = $default_vals->{'csmdata'};
  }
  my $datdir = "$rootdir";
  my $waccm_datdir = "$rootdir";

  # Plant function types.
  unless (defined($NLref->{'fpftcon'})) {
      $NLref->{'fpftcon'} = namelist::quote_string("$datdir/$default_vals->{'fpftcon'}");
  }
  $self->checkinputfile('fpftcon') if $optsref->{'test'};

  # Initial conditions
  unless ( defined($NLref->{'finidat'}) ) {
      if ( defined($default_vals->{'finidat'}) ) {
	  $NLref->{'finidat'} = namelist::quote_string("$waccm_datdir/$default_vals->{'finidat'}");
      }
  }
  if ( defined $NLref->{'finidat'} and $optsref->{'test'} ) {
      $self->checkinputfile('finidat');
  }

  # Surface datasets
  # code will use fsurdat if defined; if set to blank, code will build from raw 
  # files; if not defined and no raw file specifiers are defined, will first 
  # look for default fsurdat before building with default raw files

  # check if any raw file specifiers are defined
  if ( defined($NLref->{'mksrf_fvegtyp'}) or defined($NLref->{'mksrf_fsoitex'}) 
       or defined($NLref->{'mksrf_fsoicol'}) or defined($NLref->{'mksrf_flanwat'}) 
       or defined($NLref->{'mksrf_furban'}) or defined($NLref->{'mksrf_fglacier'}) 
       or defined($NLref->{'mksrf_flai'}) ) {
      $raw_files_spec = 1;   #true
  } else {
      $raw_files_spec = 0;   #false
  }

  if ( defined($NLref->{'fsurdat'}) ) {
      if ($NLref->{'fsurdat'} =~ /'\s*'/) {
	  delete $NLref->{'fsurdat'};   # undefine so it can build from raw files later
      } else {
	  if ( $raw_files_spec) {
	      print "clm2exp: WARNING> fsurdat already defined to be $NLref->{'fsurdat'}\n";
	      print "clm2exp: WARNING> - specified raw surface files will be ignored\n";
	      delete $NLref->{'mksrf_fvegtyp'};  #do not include in final namelist
	      delete $NLref->{'mksrf_fsoitex'};  #do not include in final namelist
	      delete $NLref->{'mksrf_fsoicol'};  #do not include in final namelist
	      delete $NLref->{'mksrf_flanwat'};  #do not include in final namelist
	      delete $NLref->{'mksrf_furban'};  #do not include in final namelist
	      delete $NLref->{'mksrf_fglacier'};  #do not include in final namelist
	      delete $NLref->{'mksrf_flai'};  #do not include in final namelist
	  }
      }
  } else {
      if (!$raw_files_spec) {
	  if ( defined($default_vals->{'fsurdat'}) ) {
             $NLref->{'fsurdat'} = namelist::quote_string("$waccm_datdir/$default_vals->{'fsurdat'}");
	  }
      }
  }

  if ( defined $NLref->{'fsurdat'} and $optsref->{'test'} ) {
      $self->checkinputfile('fsurdat');
  }

  # High resolution surface datasets
  unless ( defined($NLref->{'fsurdat'}) ) {
      unless ( defined($NLref->{'mksrf_fvegtyp'}) ) {
	  $NLref->{'mksrf_fvegtyp'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fvegtyp'}");
	  $self->checkinputfile('mksrf_fvegtyp') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_fsoitex'}) ) {
	  $NLref->{'mksrf_fsoitex'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fsoitex'}");
	  $self->checkinputfile('mksrf_fsoitex') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_fsoicol'}) ) {
	  $NLref->{'mksrf_fsoicol'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fsoicol'}");
	  $self->checkinputfile('mksrf_fsoicol') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_flanwat'}) ) {
	  $NLref->{'mksrf_flanwat'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_flanwat'}");
	  $self->checkinputfile('mksrf_flanwat') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_furban'}) ) {
	  $NLref->{'mksrf_furban'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_furban'}");
	  $self->checkinputfile('mksrf_furban') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_fglacier'}) ) {
	  $NLref->{'mksrf_fglacier'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fglacier'}");
	  $self->checkinputfile('mksrf_fglacier') if $optsref->{'test'};
      }
      unless ( defined($NLref->{'mksrf_flai'}) ) {
	  $NLref->{'mksrf_flai'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_flai'}");
	  $self->checkinputfile('mksrf_flai') if $optsref->{'test'};
      }
  }
}

#============================================================================

1   # to make use or require happy
