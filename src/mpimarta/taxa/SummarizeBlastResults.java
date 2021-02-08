package mpimarta.taxa;

import java.text.SimpleDateFormat;
import java.util.Calendar;

public class SummarizeBlastResults {

	private static SimpleDateFormat fmt;

	static {
		fmt = new SimpleDateFormat();
		fmt.setCalendar(Calendar.getInstance());
		fmt.applyPattern("MMddyy-HHmmss");
	}

	public static void main(String[] args){

		double percentile = 98;
		int percentageIdentity = 97;

		String blastResultsFileName = null; 
		String cutoffs = "1,1,1,1,1,2/3";
		String groupId = null; // TODO:: GET THE HEAD DIR

		if( args.length > 0 ){
			if( args[0].equals("--examples")){
				System.out.println("example usage (all opts):");
				System.out.println("java -Xmx128m -Xms128m -cp ~args~ marta.taxa.AssignTaxaMPI file -out= -tile= -co= -minscore -qsubmem -recall -noresolve -cutoffs= -group=");
				System.exit(0);
			}

			blastResultsFileName = args[0];

			StringBuffer options = new StringBuffer();
			for( int i = 1; i < args.length; i++ ){
				options.append(args[i]).append(" ");

			} options.append(" "); // one terminal space to match on for the last option.

			String arg = options.toString(), tmp;

			if( arg.indexOf(AssignTaxaMPI.VOTING_GROUP) != -1 ){					// a literal to use as a name.
				groupId = arg.substring(arg.indexOf(AssignTaxaMPI.VOTING_GROUP) + AssignTaxaMPI.VOTING_GROUP.length(), arg.indexOf(" ",arg.indexOf(AssignTaxaMPI.VOTING_GROUP)));
			}

			if( arg.indexOf(AssignTaxaMPI.SLIPSCORE_LIMIT) != -1){						// 
				tmp = arg.substring(arg.indexOf(AssignTaxaMPI.SLIPSCORE_LIMIT) + AssignTaxaMPI.SLIPSCORE_LIMIT.length(), arg.indexOf(" ",arg.indexOf(AssignTaxaMPI.SLIPSCORE_LIMIT)));
				percentile = Double.parseDouble(tmp);
				if( percentile < 1 ){
					percentile = percentile * 100; 
				}
			}

			if( arg.indexOf(AssignTaxaMPI.TAXONOMIC_RANK_CUTOFFS) != -1 ){		// cutoffs at each taxonomic rank.
				cutoffs = arg.substring(arg.indexOf(AssignTaxaMPI.TAXONOMIC_RANK_CUTOFFS) + AssignTaxaMPI.TAXONOMIC_RANK_CUTOFFS.length(), arg.indexOf(" ",arg.indexOf(AssignTaxaMPI.TAXONOMIC_RANK_CUTOFFS)));
			}

			if( arg.indexOf(AssignTaxaMPI.PERCENT_IDENTITY) != -1 ){							// 
				percentageIdentity = Integer.parseInt( 
						arg.substring(arg.indexOf(AssignTaxaMPI.PERCENT_IDENTITY) + AssignTaxaMPI.PERCENT_IDENTITY.length(), arg.indexOf(" ",arg.indexOf(AssignTaxaMPI.PERCENT_IDENTITY))));
			}

			// we should throw an error if the user specifies a flag/argument that we do not support.
			TaxaLogger.logInfo("Using source file: " + blastResultsFileName, null);

			try{
				TallyVotes tally = new TallyVotes();
				tally.vote( blastResultsFileName, groupId, percentile, cutoffs, percentageIdentity );
				
			} catch( Exception ex ){
				ex.printStackTrace();
				ex.printStackTrace(System.err);
				TaxaLogger.logError("Ex message:  " + ex.getMessage(), null);

			} finally {
				TaxaLogger.logInfo("The program has finished running.", null);
				System.out.println("The program has finished running.\n");
			}

		} else {
			System.err.println("Usage: marta.taxa.AssignTaxaMPI fastaFile [optional args, e.g. tophits or different outputdir]");
			System.exit(1);
		}
	}
}

