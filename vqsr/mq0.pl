#!/usr/bin/perl

while (<STDIN>)
{
    if (/^\#/)
    {
	print;
	next;
    }

    chomp;
    die $_ unless /^(.*);MQ0=(\d+)(.*)$/;
    
    print "$1;MQ0=";
    print $2+rand();
    print "$3\n";
    $max=$2 if $max< $2;
}

die unless $max;