#
#	camexp.pm			Erik Kluzek
#
#	Perl module to create a namelist for CAM.
#
#	Description of methods:
#
#	new ----------------------  Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package camexp;

use atmlndnl;
@camexp::ISA = qw(atmlndnl  namelist);

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

  my $self = $class->SUPER::new( "camexp", \%main::CAMEXP, $interactive, $file,
                                 "DefaultCAMEXPNamelist.xml", $default_vals,
                                 $CAM_config, $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  $self->{DYNAMICS}   = $CAM_config->cfg("DYNAMICS");  # model dynamics
  $self->{PERGRO}     = $CAM_config->cfg("PERGRO");    # PERGRO error-growth option
  $self->{PLEV}       = $CAM_config->cfg("PLEV");      # Ocean/sea-ice directory to use (dom/som)
  $self->{PHYSICS}    = $CAM_config->cfg("PHYSICS");   # physics directory to use
  $self->{CHEMISTRY}  = $CAM_config->cfg("CHEMISTRY"); # chemistry to use
  $self->{RESOLUTION} = $CAM_config->cfg("RESOLUTION");# horizontal resolution

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the CAM namelist variables.

  my ($self, %settings) = @_;

  my $class = ref($self);
  my $nm = "$class\:\:set_default_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $default_vals = $self->{'default_vals'};
  my $opt;

  # Get the default values from the XML file
  $self->get_default_values( %settings );

  my $waccm_chem = $self->{CHEMISTRY} eq 'waccm_mozart' ||
                   $self->{CHEMISTRY} eq 'waccm_ghg';


  unless (defined($NLref->{'nsrest'})) {
      if ( defined($default_vals->{'nsrest'}) ) {
	  $NLref->{'nsrest'} = $default_vals->{'nsrest'};
      } else {
	  die "ERROR:  Cannot determine value for nsrest\n";
      }
  }

  # Case name
  if (defined($optsref->{'case'})) {
      $opt = $optsref->{'case'};
  } elsif (defined($NLref->{'caseid'})) {
      $opt = $NLref->{'caseid'};
  } else {
      $opt = $default_vals->{'caseid'};
  }
  my $case = $opt;
  if ( $case !~ /([^\/]{1,32})([^\/]*)/ ) {
    die "$nm Bad casename: $case\n";
  }
  # If casename too long
  if ( $2 ne "" ) {
    $case = $1;
    print "\n\nWARNING:: Truncating caseid in namelist from $opt to $case\n\n" if ($self->{'printlev'});
  }
  $NLref->{'caseid'} = namelist::quote_string($case);

  # Length of simulation
  unless ( defined($NLref->{'nelapse'}) or defined($NLref->{'nestep'}) or
	   defined($NLref->{'stop_ymd'})) {
    $NLref->{'nelapse'} = $default_vals->{'nelapse'};
  }

  # Orbit
  unless ( defined($NLref->{'obliq'}) and defined($NLref->{'eccen'}) and
	   defined($NLref->{'mvelp'}) ) {
      unless ( defined($NLref->{'iyear_ad'}) ) {
	  $NLref->{'iyear_ad'} = $default_vals->{'iyear_ad'};
      }
  }

  # Timestep size
  unless (defined($NLref->{'dtime'})) {
      if ( defined($default_vals->{'dtime'}) ) {
	  $NLref->{'dtime'} = $default_vals->{'dtime'};
      }
  }

  # PERGRO settings
  if ($self->{PERGRO} eq 'TRUE') {
      unless (defined($NLref->{'prognostic_icesnow'})) {
	  $NLref->{'prognostic_icesnow'} = ".false.";
      }
  }

  # Print conservation errors
  # Turn this off for PERGRO runs or for WACCM runs
  if ($self->{PERGRO} eq 'TRUE'  or  $waccm_chem) {
      unless (defined($NLref->{'print_energy_errors'})) {
	  $NLref->{'print_energy_errors'} = ".false.";
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
  my $am_datdir = "$rootdir";

  # Initial conditions
  if ( defined($NLref->{'ncdata'}) ) {
      $opt = $NLref->{'ncdata'};
  } else {
      $opt = "$am_datdir/$default_vals->{'ncdata'}";
  }
  $NLref->{'ncdata'} = namelist::quote_string($opt);
  $self->checkinputfile('ncdata') if $optsref->{'test'};

  # Topography
  undef $opt;
  if ( defined($NLref->{'bnd_topo'}) ) {
      $opt = $NLref->{'bnd_topo'};
  } elsif (defined($default_vals->{'bnd_topo'})) {
      $opt = "$am_datdir/$default_vals->{'bnd_topo'}";
  }
  if (defined($opt)) {
      $NLref->{'bnd_topo'} = namelist::quote_string($opt);
      $self->checkinputfile('bnd_topo') if $optsref->{'test'};
  }

  # Absorptivity and emissivity data
  if ( defined($NLref->{'absems_data'}) ) {
      $opt = $NLref->{'absems_data'};
  } else {
      $opt = "$am_datdir/$default_vals->{'absems_data'}";
  }
  $NLref->{'absems_data'} = namelist::quote_string($opt);
  $self->checkinputfile('absems_data') if $optsref->{'test'};

  # SST dataset
  if ( defined($NLref->{'bndtvs'}) ) {
      $opt = $NLref->{'bndtvs'};
  } else {
      $opt = "$am_datdir/$default_vals->{'bndtvs'}";
  }
  $NLref->{'bndtvs'} = namelist::quote_string($opt);
  $self->checkinputfile('bndtvs') if $optsref->{'test'};

  # Ozone dataset
  if ( defined($NLref->{'bndtvo'}) ) {
      $opt = $NLref->{'bndtvo'};
  } else {
      $opt = "$am_datdir/$default_vals->{'bndtvo'}";
  }
  $NLref->{'bndtvo'} = namelist::quote_string($opt);
  $self->checkinputfile('bndtvo') if $optsref->{'test'};

  # Aerosol Mass climatology dataset
  if ( defined($NLref->{'bndtvaer'}) ) {
      $opt = $NLref->{'bndtvaer'};
  } else {
      $opt = "$am_datdir/$default_vals->{'bndtvaer'}";
  }
  $NLref->{'bndtvaer'} = namelist::quote_string($opt);
  $self->checkinputfile('bndtvaer') if $optsref->{'test'};

  # Aerosol optics dataset
  if ( defined($NLref->{'aeroptics'}) ) {
      $opt = $NLref->{'aeroptics'};
  } else {
      $opt = "$am_datdir/$default_vals->{'aeroptics'}";
  }
  $NLref->{'aeroptics'} = namelist::quote_string($opt);
  $self->checkinputfile('aeroptics') if $optsref->{'test'};

  # Time-variant solar constant boundary dataset
  if ( $NLref->{'scenario_scon'} =~ /RAMPED/ ) {
      if ( defined($NLref->{'bndtvscon'}) ) {
	  $opt = $NLref->{'bndtvscon'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvscon'}";
      }
      $NLref->{'bndtvscon'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvscon') if $optsref->{'test'};
  }

  # Time-variant greenhouse gas surface value boundary dataset
  if ( $NLref->{'scenario_ghg'} =~ /RAMPED/ ) {
      if ( defined($NLref->{'bndtvghg'}) ) {
	  $opt = $NLref->{'bndtvghg'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvghg'}";
      }
      $NLref->{'bndtvghg'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvghg') if $optsref->{'test'};
  }

  # Volcanic Aerosol Mass climatology dataset
  if ( $NLref->{'strat_volcanic'} =~ /\.true\./i ) {
      if ( defined($NLref->{'bndtvvolc'}) ) {
	  $opt = $NLref->{'bndtvvolc'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvvolc'}";
      }
      $NLref->{'bndtvvolc'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvvolc') if $optsref->{'test'};
  }

  # Historic carbon scaling factors.
  if ( $NLref->{'scenario_carbon_scale'} =~ /RAMPED/ ) {
#      if ( defined($NLref->{'bndtvcarbonscale'}) ) {
#	  $opt = $NLref->{'bndtvcarbonscale'};
#      } else {
#	  $opt = "$am_datdir/$default_vals->{'bndtvcarbonscale'}";
#      }
#      $NLref->{'bndtvcarbonscale'} = namelist::quote_string($opt);
#      $self->checkinputfile('bndtvcarbonscale') if $optsref->{'test'};
  }

  # carbon emissions dataset for prognostic carbon aerosols
  if ( $NLref->{'aero_carbon'} =~ /\.true\./i ) {
#      if ( defined($NLref->{'co_emis'}) ) {
#	  $opt = $NLref->{'co_emis'};
#      } else {
#	  $opt = "$am_datdir/$default_vals->{'co_emis'}";
#      }
#      $NLref->{'co_emis'} = namelist::quote_string($opt);
#      $self->checkinputfile('co_emis') if $optsref->{'test'};
  }

  # Prognostic sulfur cycle.
  if ( $NLref->{'prognostic_sulfur'} =~ /direct/ or $NLref->{'prognostic_sulfur'} =~ /passive/ ) {

      # DMS emissions dataset
      if ( defined($NLref->{'bndtvdms'}) ) {
	  $opt = $NLref->{'bndtvdms'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvdms'}";
      }
      $NLref->{'bndtvdms'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvdms') if $optsref->{'test'};

      # oxidant dataset
      if ( defined($NLref->{'bndtvoxid'}) ) {
	  $opt = $NLref->{'bndtvoxid'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvoxid'}";
      }
      $NLref->{'bndtvoxid'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvoxid') if $optsref->{'test'};

      # SOx emissions dataset
      if ( defined($NLref->{'bndtvsox'}) ) {
	  $opt = $NLref->{'bndtvsox'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvsox'}";
      }
      $NLref->{'bndtvsox'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvsox') if $optsref->{'test'};
  }

  # soil erodibility dataset
#  if ( defined($NLref->{'soil_erod'}) ) {
#      $opt = $NLref->{'soil_erod'};
#  } else {
#      $opt = "$am_datdir/$default_vals->{'soil_erod'}";
#  }
#  $NLref->{'soil_erod'} = namelist::quote_string($opt);
#  $self->checkinputfile('soil_erod') if $optsref->{'test'};

  # simulated ISCCP dataset
  if ( $NLref->{'doisccp'} =~ /\.true\./i ) {
      if ( defined($NLref->{'isccpdata'}) ) {
	  $opt = $NLref->{'isccpdata'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'isccpdata'}";
      }
      $NLref->{'isccpdata'} = namelist::quote_string($opt);
      $self->checkinputfile('isccpdata') if $optsref->{'test'};
  }

  # Greenhouse gas production/loss rates
  if ( $NLref->{'trace_gas'} =~ /t/i  or $self->{CHEMISTRY} eq 'waccm_ghg' ) {
      if ( defined($NLref->{'bndtvg'}) ) {
	  $opt = $NLref->{'bndtvg'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvg'}";
      }
      $NLref->{'bndtvg'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvg') if $optsref->{'test'};
  }

  # WACCM options
  if ($waccm_chem) {

      # O2,O1,N2, CO2 Constituents for non-LTE calculations and heating rates below 200 nm
      if ( defined($NLref->{'cftgcm'}) ) {
	  $opt = $NLref->{'cftgcm'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'cftgcm'}";
      }
      $NLref->{'cftgcm'} = namelist::quote_string($opt);
      $self->checkinputfile('cftgcm') if $optsref->{'test'};

  }

  # WACCM_GHG options
  if ( $self->{CHEMISTRY} eq 'waccm_ghg' ) {
      if ( defined($NLref->{'h2orates'}) ) {
	  $opt = $NLref->{'h2orates'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'h2orates'}";
      }
      $NLref->{'h2orates'} = namelist::quote_string($opt);
      $self->checkinputfile('h2orates') if $optsref->{'test'};
  }

  # WACCM_MOZART options
  if ( $self->{CHEMISTRY} eq 'waccm_mozart' ) {

      if ( defined($NLref->{'lbc_file'}) ) {
	  $opt = $NLref->{'lbc_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'lbc_file'}";
      }
      $NLref->{'lbc_file'} = namelist::quote_string($opt);
      $self->checkinputfile('lbc_file') if $optsref->{'test'};

      if ( defined($NLref->{'chem_config'}) ) {
	  $opt = $NLref->{'chem_config'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'chem_config'}";
      }
      $NLref->{'chem_config'} = namelist::quote_string($opt);
      $self->checkinputfile('chem_config') if $optsref->{'test'};

      if ( defined($NLref->{'airpl_emis_file'}) ) {
	  $opt = $NLref->{'airpl_emis_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'airpl_emis_file'}";
      }
      $NLref->{'airpl_emis_file'} = namelist::quote_string($opt);
      $self->checkinputfile('airpl_emis_file') if $optsref->{'test'};

      if ( defined($NLref->{'nox_emis_file'}) ) {
	  $opt = $NLref->{'nox_emis_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'nox_emis_file'}";
      }
      $NLref->{'nox_emis_file'} = namelist::quote_string($opt);
      $self->checkinputfile('nox_emis_file') if $optsref->{'test'};

#      if ( defined($NLref->{'co_emis_file'}) ) {
#	  $opt = $NLref->{'co_emis_file'};
#      } else {
#	  $opt = "$am_datdir/$default_vals->{'co_emis_file'}";
#      }
#      $NLref->{'co_emis_file'} = namelist::quote_string($opt);
#      $self->checkinputfile('co_emis_file') if $optsref->{'test'};

      if ( defined($NLref->{'ch2o_emis_file'}) ) {
	  $opt = $NLref->{'ch2o_emis_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'ch2o_emis_file'}";
      }
      $NLref->{'ch2o_emis_file'} = namelist::quote_string($opt);
      $self->checkinputfile('ch2o_emis_file') if $optsref->{'test'};

      if ( defined($NLref->{'sad_file'}) ) {
	  $opt = $NLref->{'sad_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'sad_file'}";
      }
      $NLref->{'sad_file'} = namelist::quote_string($opt);
      $self->checkinputfile('sad_file') if $optsref->{'test'};

      if ( defined($NLref->{'sulf_file'}) ) {
	  $opt = $NLref->{'sulf_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'sulf_file'}";
      }
      $NLref->{'sulf_file'} = namelist::quote_string($opt);
      $self->checkinputfile('sulf_file') if $optsref->{'test'};

      if ( defined($NLref->{'depvel_file'}) ) {
	  $opt = $NLref->{'depvel_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'depvel_file'}";
      }
      $NLref->{'depvel_file'} = namelist::quote_string($opt);
      $self->checkinputfile('depvel_file') if $optsref->{'test'};

      if ( defined($NLref->{'n2d_file'}) ) {
	  $opt = $NLref->{'n2d_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'n2d_file'}";
      }
      $NLref->{'n2d_file'} = namelist::quote_string($opt);
      $self->checkinputfile('n2d_file') if $optsref->{'test'};

      if ( defined($NLref->{'xs_coef_file'}) ) {
	  $opt = $NLref->{'xs_coef_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'xs_coef_file'}";
      }
      $NLref->{'xs_coef_file'} = namelist::quote_string($opt);
      $self->checkinputfile('xs_coef_file') if $optsref->{'test'};

      if ( defined($NLref->{'xs_short_file'}) ) {
	  $opt = $NLref->{'xs_short_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'xs_short_file'}";
      }
      $NLref->{'xs_short_file'} = namelist::quote_string($opt);
      $self->checkinputfile('xs_short_file') if $optsref->{'test'};

      if ( defined($NLref->{'xs_long_file'}) ) {
	  $opt = $NLref->{'xs_long_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'xs_long_file'}";
      }
      $NLref->{'xs_long_file'} = namelist::quote_string($opt);
      $self->checkinputfile('xs_long_file') if $optsref->{'test'};

      if ( defined($NLref->{'rsf_file'}) ) {
	  $opt = $NLref->{'rsf_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'rsf_file'}";
      }
      $NLref->{'rsf_file'} = namelist::quote_string($opt);
      $self->checkinputfile('rsf_file') if $optsref->{'test'};

      if ( defined($NLref->{'tgcm_ubc_file'}) ) {
	  $opt = $NLref->{'tgcm_ubc_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'tgcm_ubc_file'}";
      }
      $NLref->{'tgcm_ubc_file'} = namelist::quote_string($opt);
      $self->checkinputfile('tgcm_ubc_file') if $optsref->{'test'};

      if ( defined($NLref->{'snoe_ubc_file'}) ) {
	  $opt = $NLref->{'snoe_ubc_file'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'snoe_ubc_file'}";
      }
      $NLref->{'snoe_ubc_file'} = namelist::quote_string($opt);
      $self->checkinputfile('snoe_ubc_file') if $optsref->{'test'};

  }


}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
