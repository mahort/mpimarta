/*
 * by Matt Horton
 * 
 * AssignTaxa is the 'main' application of the software bundle 'marta'. Through AssignTaxa one 
 * taxonomically classifies each sequence in a tab delimited file (containing an id and sequence to BLAST.
 * 
 * AssignTaxa creates 'jobs' using the Quartz utility. Each job is responsible for its own sequence.
 * To change voting parameters for a set of sequences that have already been investigated, see the Revote class.
 */

package mpimarta.taxa;

import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import org.apache.log4j.Logger;

public class AssignTaxaMPI {

	static Logger log4j = Logger.getLogger("mpimarta.taxa");

	protected static int[] dim( String filename, char delim ) throws IOException {

		Reader reader = new InputStreamReader( new FileInputStream( filename ));        

        char[] buffer = new char[4096];
        int line = 0, field = 0;

        for( int read = reader.read( buffer ); read >= 0; read = reader.read( buffer )){

            for( int i = 0; i < read ; i++){
                if(( line == 0 ) && ( buffer[i] == delim )){
                	field++;

                } else if( buffer[i] == '\n' ){
                	line++;
                }
            }

        } field++; // get the last field, since counting the delimiters only gives us the (n-1) separators.

        reader.close();
		System.gc();

        return new int[]{line, field};
	}

	private static final String BLAST_ONLY = "-nur";
	private static final String PROJECT_NAME = "-proj=";
	private static final String WALL_TIME = "-wall=";
	private static final String NCBI_DB = "-db=";
	private static final String NUMBER_OF_CORES = "-co=";
	private static final String BLAST_WORDSIZE = "-ws=";

	protected static final String MINSCORE = "-minscore";
	protected static final String NUMBER_OF_HITS_TO_CONSIDER = "-top=";
	protected static final String PERCENT_IDENTITY = "-per=";
	protected static final String SLIPSCORE_LIMIT = "-tile=";
	protected static final String TAXONOMIC_RANK_CUTOFFS = "-cutoffs=";
	protected static final String VOTING_GROUP = "-group=";

	public static void main(String[] args){

		boolean minscore = false, blastOnly = false;

		double percentile = 98;
		int percentageIdentity = 97;
		int tophits = 100, wordsize = 28, cores = 48;

		String fastaFileName = null; 
		String cutoffs = "1,1,1,1,1,2/3";
		String groupId = null, projectName = "cegs", wallTime="08:00:00", ncbiDatabase = "nt", outputdir = "~/"; // TODO:: GET THE HEAD DIR

		if( args.length > 0 ){
			if( args[0].equals("--examples")){
				System.out.println("example usage (all opts):");
				System.out.println("java -Xmx128m -Xms128m -cp ~args~ marta.taxa.AssignTaxaMPI file -out= -tile= -co= -minscore -qsubmem -recall -noresolve -cutoffs= -group=");
				System.exit(0);
			}

			fastaFileName = args[0];
			outputdir = fastaFileName.substring(0, fastaFileName.lastIndexOf("/"));			

			StringBuffer options = new StringBuffer();
			for( int i = 1; i < args.length; i++ ){
				options.append(args[i]).append(" ");

			} options.append(" "); // one terminal space to match on for the last option.

			String arg = options.toString(), tmp;
//
// booleans
			minscore = ( arg.indexOf(MINSCORE) != -1 );				// 
			blastOnly = ( arg.indexOf( BLAST_ONLY ) != -1 );		// don't vote (e.g. mysql is installed elsewhere).

			if( arg.indexOf(NUMBER_OF_HITS_TO_CONSIDER) != -1){		// the number of hits to consider; I typically use 100; 50 is probably fine.
				tmp = arg.substring(arg.indexOf(NUMBER_OF_HITS_TO_CONSIDER) + NUMBER_OF_HITS_TO_CONSIDER.length(), arg.indexOf(" ",arg.indexOf(NUMBER_OF_HITS_TO_CONSIDER)));
				tophits = Integer.parseInt(tmp);
			}

			if( arg.indexOf(NCBI_DB) != -1){						// e.g. nt
				ncbiDatabase = arg.substring(arg.indexOf(NCBI_DB) + NCBI_DB.length(), arg.indexOf(" ",arg.indexOf(NCBI_DB)));
			}

			if( arg.indexOf(NUMBER_OF_CORES) != -1){				// the number of cores to use for mpi (for mpiblast: # of db fragments + 2)
				tmp = arg.substring(arg.indexOf(NUMBER_OF_CORES) + NUMBER_OF_CORES.length(), arg.indexOf(" ", arg.indexOf(NUMBER_OF_CORES)));
				cores = Integer.parseInt(tmp);
			}

			if( arg.indexOf(WALL_TIME) != -1){				// the number of cores to use for mpi (for mpiblast: # of db fragments + 2)
				wallTime = arg.substring(arg.indexOf(WALL_TIME) + WALL_TIME.length(), arg.indexOf(" ", arg.indexOf(WALL_TIME)));
			}

			if( arg.indexOf(SLIPSCORE_LIMIT) != -1){						//  
				tmp = arg.substring(arg.indexOf(SLIPSCORE_LIMIT) + SLIPSCORE_LIMIT.length(), arg.indexOf(" ",arg.indexOf(SLIPSCORE_LIMIT)));
				percentile = Double.parseDouble(tmp);
			}

			if( arg.indexOf(VOTING_GROUP) != -1 ){					// a literal to use as a name.
				groupId = arg.substring(arg.indexOf(VOTING_GROUP) + VOTING_GROUP.length(), arg.indexOf(" ",arg.indexOf(VOTING_GROUP)));
			}

			if( arg.indexOf(BLAST_WORDSIZE) != -1){					// wordsize, defaults to 28 like megablast.
				tmp = arg.substring(arg.indexOf(BLAST_WORDSIZE) + BLAST_WORDSIZE.length(), arg.indexOf(" ",arg.indexOf(BLAST_WORDSIZE)));
				wordsize = Integer.parseInt(tmp);
			}

			if( arg.indexOf(TAXONOMIC_RANK_CUTOFFS) != -1 ){		// cutoffs at each taxonomic rank.
				cutoffs = arg.substring(arg.indexOf(TAXONOMIC_RANK_CUTOFFS) + TAXONOMIC_RANK_CUTOFFS.length(), arg.indexOf(" ",arg.indexOf(TAXONOMIC_RANK_CUTOFFS)));
			}

			if( arg.indexOf(PERCENT_IDENTITY) != -1 ){							// 
				percentageIdentity = Integer.parseInt( 
						arg.substring(arg.indexOf(PERCENT_IDENTITY) + PERCENT_IDENTITY.length(), arg.indexOf(" ",arg.indexOf(PERCENT_IDENTITY))));
			}

			if( arg.indexOf(PROJECT_NAME) != -1 ){							// 
				projectName = arg.substring(arg.indexOf(PROJECT_NAME) + PROJECT_NAME.length(), arg.indexOf(" ",arg.indexOf(PROJECT_NAME)));
			}

			// we should throw an error if the user specifies a flag/argument that we do not support.

		} else {
			System.err.println("Usage: marta.taxa.AssignTaxaMPI fastaFile [optional args, e.g. tophits or different outputdir]");
			System.exit(1);
		}

		log4j.info("Using source file: " + fastaFileName );
		log4j.info("Using output directory: " + outputdir );
		log4j.info("Running AssignTaxaMPI with minscore option set to: " + minscore );

// parse the file, assuming that it is in the following format:
// first column: sequence-name
// second column: sequence

		try{

// split the file up and recall with qsub.
			String scriptName = writeQsubScript( projectName, wallTime, fastaFileName, ncbiDatabase, minscore, cores, 
					groupId, cutoffs, percentile, percentageIdentity, wordsize, tophits, blastOnly );

			if( null == scriptName ){
				System.err.println("Unable to partition the blast file into smaller partitions for cluster-sequencing");
				System.exit(1);
			}

			boolean success = fireQsubScript( scriptName );
			if( !success ){
				System.err.println("Unable to queue the blast jobs.");
				System.exit(1);
			}

			return;

		} catch( Exception ex ){
			log4j.error( "Error occurred at "+ getTime());
			ex.printStackTrace();
			ex.printStackTrace(System.err);
			log4j.error("Ex message:  " + ex.getMessage());

		} finally {}

		return;
	}

	private static String writeQsubScript( String projectName, String wallTime, String filename, String ncbiDatabase, boolean minscore, int numberOfCores, 
			String groupId, String cutoffs, double slippageTolerance, double percentIdentity, int wordSize, int topHits, boolean skipVote ){

		StringBuffer qsubCmds = new StringBuffer();
		qsubCmds.append( "#!/bin/bash\n" );
		qsubCmds.append( "#PBS -S /bin/bash\n" );
		qsubCmds.append( "#PBS -N marta\n" );
		qsubCmds.append( "#PBS -j oe\n" );
		qsubCmds.append( "#PBS -P " ).append( projectName ).append( "\n");
		qsubCmds.append( "#PBS -l walltime=" ).append( wallTime ).append( "\n" );
		qsubCmds.append( "#PBS -l ncpus=48,mem=2gb\n" );
		qsubCmds.append( "#PBS -o marta.mpiblast.log\n" );

		qsubCmds.append( "\n" );
		qsubCmds.append( "cd ${PBS_O_WORKDIR}\n" );
		qsubCmds.append( "blast_file=" ).append( filename ).append( '\n' );
		qsubCmds.append( "blast_output=" ).append( filename ).append("_out\n" );
		qsubCmds.append( "mpiexec -n " ).append( numberOfCores ).append( 
						" mpiblast -m 8 -d " ).append( ncbiDatabase ).append(
						" -I T -v " ).append( topHits ).append( " -b " ).append( topHits ).append( 
						" -W " ).append( wordSize ).append( " -p blastn -i ${blast_file} -o ${blast_output_file}\n" );

		if( skipVote ){
			qsubCmds.append( "#" );
		}

		qsubCmds.append("java -Xmx128m -Xms128m -cp `cat cp.txt` mpimarta.taxa.SummarizeBlastResults ").append( 
				filename ).append( "_out" ).append(( minscore ? " " + MINSCORE : "" )).append( 
				" " + VOTING_GROUP).append( groupId ).append(
				" " + SLIPSCORE_LIMIT ).append( slippageTolerance ).append(
				( cutoffs == null ? "" : " " + TAXONOMIC_RANK_CUTOFFS.concat( cutoffs ))).append( 
				" " + PERCENT_IDENTITY ).append( percentIdentity ).toString();

		String qsubFile = new StringBuffer( filename ).append( ".sh" ).toString();

		try{

	        BufferedWriter writer = new BufferedWriter( new OutputStreamWriter( new FileOutputStream( qsubFile, false), "UTF-8"));

	        writer.write(qsubCmds.toString());
	        writer.write('\n');
	        writer.flush();
	        writer.close();

		} catch( IOException ioEx ){
			log4j.error("IOException during java blasting for file: " + filename );
			log4j.error(ioEx.getMessage());
			ioEx.printStackTrace(System.err);
			System.exit(1);
			return null;
		}

		return qsubFile;
	}

	// 
	private static boolean fireQsubScript( String scriptName ){

		try{

			String[] cmds = {
					"sh",
					"-c", // how do I delay this... is it `-a yadda`???
					new StringBuffer( "qsub " ).append( scriptName ).toString()
			};

			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec(cmds);

			InputStream err = proc.getErrorStream();
	        InputStreamReader reader = new InputStreamReader(err);
	        BufferedReader breader = new BufferedReader(reader);
	        String line; boolean error = false;
	        while(( line = breader.readLine()) != null){
	        	error = true;
	        	log4j.error(line);
	        }

	        if( error ){
	        	return false;
	        }

	        InputStream out = proc.getInputStream();
	        reader = new InputStreamReader(out);
	        breader = new BufferedReader(reader);
	        while(( line = breader.readLine()) != null){
	        	log4j.info(line);
	        }

	        try{
				proc.waitFor();
				log4j.info("TODO: update/ensure the *.sh is actually in the queue. proc.waitFor() won't do anything here!: " + scriptName );

			}  catch ( InterruptedException e ) { 
				log4j.error("InterruptedException in launchBlastJobs()!");
				e.printStackTrace(System.err);
				return false;
			}

			return true;

		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
		}

		return false;
	}

	private static String getTime(){
		SimpleDateFormat fmt = new SimpleDateFormat();
		fmt.setCalendar(Calendar.getInstance());
		fmt.applyPattern("MMddyy-HHmmss");
		return fmt.format( fmt.getCalendar().getTime());
	}
}