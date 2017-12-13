#!/usr/bin/python

## 2008/04/10
## Author: Yih-En Andrew Ban

## Sort lines in a text file by multiple field comparison in order of field importance (highest -> lowest).  Handles numeric data only.

import sys
from optparse import OptionParser

########
# main #
########
def main( options, args ):
	# parse fields, no re-indexing from zero because we need negative
	# numbers (i.e. -0 doesn't exist)
	fields = options.fields.strip().split( ',' )
	fields = [ int( i ) for i in fields ]

	# maps tuple key -> line
	data = {}

	# open input
	if len( args ) > 0:
		f = open( args[ 0 ] )
	else:
		f = sys.stdin

	# skip and print any header lines
	for i in range( options.skip_header ):
		line = f.readline()
		if options.print_header:
			sys.stdout.write( line )

	# run through lines and store
	for line in f:
		entry = line.strip().split()

		if len( entry ) > 0:
			if options.skip_char == None or entry[ 0 ][ 0 ] != options.skip_char:
				# store tuple key -> line
				key = []
				try:
					for i in fields:
						if i > 0:
							key.append( float( entry[ i - 1 ] ) )
						else:
							key.append( -float( entry[ -i - 1 ] ) )
				except ValueError:
					sys.stderr.write( "Error: Cannot convert field to number!  Did you forgot to skip header lines?  Line was:\n%s" % ( line ) )
					sys.exit( 1 )

				data[ tuple( key ) ] = line

	# close input only if not stdin
	if len( args ) > 0:
		f.close()

	# sort tuple keys
	keys = data.keys()
	keys.sort()
	if options.reverse_sort:
		keys.reverse()

	# output
	for k in keys:
		sys.stdout.write( data[ k ] )


##############
# invocation #
##############
if __name__ == '__main__':
	usage = "usage: %prog [OPTIONS] input"
	parser = OptionParser( usage, version=0.001 )
	parser.description = "Sort lines in a text file by multiple field comparison in order of field importance (highest -> lowest).  Handles numeric data only."
	parser.add_option( "-s", "--skip_header", action="store", type="int", dest="skip_header", default=0, help="Skip N lines from beginning of file." )
	parser.add_option( "-p", "--print_header", action="store_true", dest="print_header", default=False, help="In sorted output, print the header lines that were skipped using -s/--skip_header." )
	parser.add_option( "-c", "--skip_char", action="store", type="string", dest="skip_char", default=None, help="Skip lines beginning with given character." )
	parser.add_option( "-r", "--reverse", action="store_true", dest="reverse_sort", default=False, help="Reverse the entire sort." )
	parser.add_option( "-f", "--fields", action="store", dest="fields", default=None, help="List of fields to use in order of highest to lowest importance (counting starts from 1, not 0).  Separate numbers by commas, e.g. 1,9,-4.  Negative number indicates reversal of sort for that particular field." )

	( options, args ) = parser.parse_args()

	if len( args ) == 0 and options.fields == None:
		parser.print_help()
		sys.exit( 1 )
	if options.fields == None:
		sys.stderr.write( "Error: Missing fields!\n")
		sys.exit( 1 ) 
	else:
		main( options, args )

