#Function: runs a list of commands in parallel, always occupying N processes
#input: 1. file name with the list of processes
#		2. number of processes to run in parallel
use strict; 
use lib "/home/alu/binknis/Perl_DNAE";
use lib "$ENV{HOME}/Perl_DNAE"; 
use getAll; 
use POSIX ":sys_wait_h";

(my $cmd_list_file, my $num_processes ) = @ARGV; 
my @tasks = @{ getAll::lines($cmd_list_file) }; 
my %kids;

{
  while (@tasks and keys %kids < $num_processes) {
    $kids{fork_a_task(shift @tasks)} = "active";
  }
  {
    my $pid = waitpid(-1, 0); #wait for any child process
    if ($pid == -1) { #no child process
      %kids = ();
    } else { #a child finished (any child)
      delete $kids{$pid};
    }
  }
  redo if @tasks or %kids;
}

sub fork_a_task {
  my $cmd = shift;
  my $pid = fork;
  return $pid if $pid; ### parent - return PID ###
  #### child ###
  unless (defined $pid) { #fork failed
    warn "cannot fork: $!";
    return 0;
  }
  #What to do in child - START
  $cmd =~ s/\(/\\\(/g; $cmd =~ s/\)/\\\)/g; #replace parenthesis for compatability with UNIX syntax
  system ($cmd);
  #What to do in child - END
  exit 0;
}

### code to use to check this algorithm (just prints the commands with varying sleep time)
# my $sleep_time = $$ % 30; 
  # print $sleep_time ." for ". $$ ."\n"; 
  # sleep $sleep_time; #***
  #system ("echo $cmd"); 
