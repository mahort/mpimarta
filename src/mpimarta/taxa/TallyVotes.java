package mpimarta.taxa;

import java.io.*;
import java.util.*;

public class TallyVotes {

	private HashMap<String, ArrayList<BlastCandidate>> blastResults = null;
	private HashMap<String, String> blastResultSequences = null;
	private HashMap<String, Integer> blastResultSequenceLength = null;
	private HashMap<String, Double> maxBitScores = null;

	private boolean parseMPIblastOutput( String blastResultsFileName ) throws IOException {

		String line;

		BufferedReader reader = null;
		this.blastResults = new HashMap<String, ArrayList<BlastCandidate>>();

		try{

			reader = new BufferedReader( new InputStreamReader( 
					new FileInputStream( blastResultsFileName ), "UTF-8" ));

			while(( line = reader.readLine()) != null ){

				String[] fields = line.replaceAll("  +", "\t").split("\t");
				String otuId = new StringBuffer( "otu" ).append( fields[0] ).toString();

				ArrayList<BlastCandidate> votes = null;
				if( this.blastResults.containsKey(otuId)){
					votes = this.blastResults.get(otuId);

				} else {
					votes = new ArrayList<BlastCandidate>();
				}

				BlastCandidate candidate = new BlastCandidate(otuId, fields);
				votes.add( candidate );
				this.blastResults.put(otuId, votes);
			}

			return true;
	
		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
			return false;

		} catch( Exception ex ){
			ex.printStackTrace(System.err);
			return false;

		} finally {

			if( reader != null ){
				reader.close();
			}
		}
	}

	/*
	 * we need to know the length of each target-sequence.
	 */
	private boolean filterResultsOnTargetLength( String fastaFileName, double overlap ){

		BufferedReader reader = null;

		try{
			reader = new BufferedReader( new InputStreamReader( new FileInputStream( fastaFileName ), "UTF-8" ));

			String line, otuId = new String();
			StringBuilder blastSequence = new StringBuilder();

			this.blastResultSequences = new HashMap<String, String>();
			this.blastResultSequenceLength = new HashMap<String, Integer>();

			while(( line = reader.readLine()) != null ){

				// this is a fasta file, so we're either reading an id (prefixed with >) or a sequence
				if( line.startsWith(">")){
					otuId = new StringBuffer( "otu" ).append( line.replaceFirst("^>", "").split(" ")[0] ).toString();
					blastSequence = new StringBuilder();

				} else {
					blastSequence.append( line );
				}

				this.blastResultSequences.put( otuId, blastSequence.toString());
				this.blastResultSequenceLength.put( otuId, blastSequence.toString().length());
			}

			Iterator<Map.Entry<String, ArrayList<BlastCandidate>>> entries = this.blastResults.entrySet().iterator();
			ArrayList<BlastCandidate> candidates;
			int necessaryCoverage;

			while( entries.hasNext()){
				Map.Entry<String, ArrayList<BlastCandidate>> entry = entries.next();
				otuId = entry.getKey();
				candidates = entry.getValue();

				// the original sequence length * the overlap that we expect demarcates the NECESSARY overlap (coverage) we expect.
				necessaryCoverage = (int)(((double)overlap) * this.blastResultSequenceLength.get( otuId ));
				// if there's a null pointer here, it means that the results file had an 'OTU-id' that wasn't in the fasta file. 
				TaxaLogger.logDebug( "Necessary coverage for otu " + otuId + " is " + necessaryCoverage, null);

				for( int i = (candidates.size() - 1 ); i >= 0; i-- ){
					BlastCandidate candidate = candidates.get(i);
					if( candidate.getLength() < necessaryCoverage ){
						TaxaLogger.logDebug( "Analyzing " + otuId + ", removing " + candidate.getId() + " (It is short).", null);
						candidates.remove( candidate );
					}
				}
			}

			return true;

		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
			TaxaLogger.logError("IOException in AssignTaxa", null);

		} catch( Exception ex ){
			ex.printStackTrace();
			ex.printStackTrace(System.err);
			TaxaLogger.logError("Ex message:  " + ex.getMessage(), null);

		} finally {}

		return false;
	}

	/* by default, require 80% overlap */
	private boolean filterResultsOnTargetLength( String fastaFileName ){
		return filterResultsOnTargetLength( fastaFileName, 0.80);
	}

	private boolean determineTopBitScores(){

		try{
			Iterator<Map.Entry<String, ArrayList<BlastCandidate>>> entries = this.blastResults.entrySet().iterator();

			String otuId;
			ArrayList<BlastCandidate> candidates; BlastCandidate candidate;

			this.maxBitScores = new HashMap<String, Double>();

			while( entries.hasNext()){

				Map.Entry<String, ArrayList<BlastCandidate>> entry = entries.next();
				otuId = entry.getKey();
				candidates = entry.getValue();

				double winningScore = -1, score = -1; // check to see if setting the default to 1 ensures that all candidates will be included (i.e. that nothing is omitted???)
				// iterate to get each blast-candidate's bitscore.
				for( int i = 0; i <= ( candidates.size() - 1); i++ ){
					candidate = candidates.get( i );
					score = candidate.getBitscore();

					if( Double.compare( score, winningScore ) > 0 ){
						winningScore = score;
					}
				}

				this.maxBitScores.put( otuId, new Double( winningScore ));
			}

			return true;

		} catch(Exception ex){
			ex.printStackTrace(System.err);
		}

		return false;
	}

	private boolean identifyTheWinners( String outputFileNamePrefix, String groupId, double slipScoreTolerance, String taxonomicCutoffs, double percentageIdentity ) throws IOException {

		BufferedWriter writer = null;

		try{

			// 0.) prepare to write the results
			writer = new BufferedWriter(
					new OutputStreamWriter( new FileOutputStream( 
					new StringBuffer( outputFileNamePrefix ).append( "_taxa.summary" ).append(
												( groupId == null ? ".txt" : "." + groupId.trim() + ".txt")).toString(), false), "UTF-8"));

	        writer.write( "fasta_id\tgi\ttaxon_id\ttaxonomic_rank\taffinity\te_value\twinning_score\tmax_score\tperc_id\tcoverage\tfull_lineage\tvotes\toverall\n");

	        // 1.) sort the collections. Print these to see if they need to be reverse-sorted (I'm sure they do).
			Iterator<Map.Entry<String, ArrayList<BlastCandidate>>> entries = this.blastResults.entrySet().iterator();
			ArrayList<BlastCandidate> candidates;
			String otuId;

			while( entries.hasNext()){
				Map.Entry<String, ArrayList<BlastCandidate>> entry = entries.next();
				otuId = entry.getKey();
				candidates = entry.getValue();

				Collections.sort(candidates);
				this.blastResults.put( otuId, candidates);
			}

			// 2.) iterate over the candidates and vote.
			int candidateCount = 0;
			entries = this.blastResults.entrySet().iterator();
			while( entries.hasNext()){

				Map.Entry<String, ArrayList<BlastCandidate>> entry = entries.next();
				otuId = entry.getKey();
				candidates = entry.getValue();

				candidateCount++;
				if( candidateCount % 500 == 0 ){
					System.out.println("On candidate: " + candidateCount + " of " + this.blastResults.size() + "\n");
				}

				///////////////////////////////////////////////////////////////////////////////////////////////
				//A.) get the taxonomic information (lineage, taxaIds) from our ncbi database.
				Hashtable<String, Taxonomy> taxonomicInfo = TaxaDatabase.getTaxonomicInformation( otuId, candidates );

				///////////////////////////////////////////////////////////////////////////////////////////////
				// summarize the tally for the valid reads...
				///////////////////////////////////////////////////////////////////////////////////////////////
				double maximumBitScore = this.maxBitScores.get(otuId).doubleValue();
				double slipScore = maximumBitScore * ( slipScoreTolerance / 100 );
				TaxaLogger.logDebug( "Using slipscore: " + slipScore, otuId );
				// when I voted by min-evalue the argument was::: double slipScore = maximumBitScore / getSlippageTolerance();

				Vector<String[]> votes = new Vector<String[]>();
				Vector<String[]> votesForARainyDay = new Vector<String[]>();
				double bitScore = -1;
				double eValue = -1;

				///////////////////////////////////////////////////////////////////////////////////////////////
				//A.1) some GIs will not have taxonomic information in the database. Remove them from the array.
				for( int i = (candidates.size() - 1 ); i >= 0; i-- ){
					BlastCandidate candidate = candidates.get(i);
					if( !taxonomicInfo.containsKey( candidate.getId())){
						TaxaLogger.logError("No taxonomic information for GI: " + candidate.getId() + " with score: " + candidate.getBitscore(), otuId);
						candidates.remove(i);
					}
				}

				// THIS FOR STATEMENT 'VOTES' AND HANDLES THE SLIP-SCORE STRATEGY (CONTINUING AMONG BITSCORE TIES)
				// AIMS:1) step down the levels (from the 'top' score to the slip-score (basement))
				// AIMS:2) record the taxonomic ranks; if the GI values aren't labeled 'uncultured' (or something similar), vote with the information we do have.
				// B. ITERATE OVER CANDIDATES (SORTED BY BITSCORE)
				for( int j = 0; j <= ( candidates.size() - 1); j++ ){

					BlastCandidate candidate_j = candidates.get(j);

					///////////////////////////////////////////////////////////////////////////////////////////////
					// B.1) Determine the gi's Taxonomic/lineage info
					///////////////////////////////////////////////////////////////////////////////////////////////
					String gi = candidate_j.getId();
					bitScore = candidate_j.getBitscore();
					eValue = candidate_j.getMinimumEvalue();

					///////////////////////////////////////////////////////////////////////////////////////////////
					// B.2) TAXONOMIC INFO FROM GI
					// warning: not-uncommon for a null pointer @ this line when there's a disconnect b/w the database and the local blast database.
					// seems to happen even if the database is BRAND SPANKING NEW FROM THE WEBSITE b/c the taxonomy info is only updated 1x a week
					// it is also possible that the GIs are available in the blast-utility supporting files, but not the taxonomy database, because
					// of errors during the curation of the databases at ncbi.
					///////////////////////////////////////////////////////////////////////////////////////////////
					// We use the GIs (from the blast results) to determine the taxon-ids (using the gi to get the taxonid and literal info).
					// note! The 'most recently analyzed GI', at this level, may not have the winning taxonomicId (bien sur).
					Taxonomy taxonomy_j = taxonomicInfo.get( gi ); 

					///////////////////////////////////////////////////////////////////////////////////////////////
					// B.3) Record the 'vote' from the gi-taxonomic information.
					// However: SKIP uncultured/unidentified critters. One clue: no genus level assignment!
					// CORRECTION: THERE ARE SOME SPP (HITHERTO UNCERTAIN) THAT DO NOT HAVE GENUS LEVEL ASSIGNMENTS, BUT
					//	 				WHICH SEEM TO BE ACCURATE, GIVEN THE CORPUS TEST.
					// QUESTION/TO TEST-CASE: ARE THERE SAMPLES W/ MISSING GENERA AND SPP NODES?
					//					THE TAXONOMY DATABASE CONTINUES TO BE CURATED.
					///////////////////////////////////////////////////////////////////////////////////////////////
					if( null == taxonomy_j.getGenus() || 
							taxonomy_j.getSpecies().toLowerCase().startsWith("Uncultur".toLowerCase()) ||
							taxonomy_j.getSpecies().toLowerCase().startsWith("Unknown".toLowerCase()) || 
							taxonomy_j.getSpecies().toLowerCase().startsWith("Unidentified".toLowerCase())){
						
						TaxaLogger.logDebug("GI: " + gi + " lacks usable taxonomic information. Score-level: " + bitScore, otuId);

						// Only record the best results...
						if( Double.compare( bitScore, maximumBitScore) == 0 ){
							votesForARainyDay.add( taxonomy_j.getRestrictedTaxonIds());
						}

					} else {
						// wtf does this mean ---> we only return "name_class='scientific name'" which should be 'having n = 1' for each TaxonId per my 'group by'
						// we vote... by TAXONID. <-- c'est vrai.
						votes.add( taxonomy_j.getTaxonIds());
					}

					// B.4) The slip-score logic:
					// B.4.1) If there are ties, continue:
					if(( j != (candidates.size() - 1 )) && ( Double.compare( bitScore, Double.valueOf( candidates.get(j+1).getBitscore())) == 0 )){
						// there are still ties to analyze.
						continue;

					} else { // B.4.2) There are NO MORE TIES to consider...

						if( votes.size() > 0 ){
							///////////////////////////////////////////////////////////////////////////////////////////////
							// We have votes and there are no more records at this 'level', so summarize what we have.
							///////////////////////////////////////////////////////////////////////////////////////////////
							break;

						} else {
							if(( j != (candidates.size() - 1)) && ( Double.compare( slipScore, candidates.get(j+1).getBitscore()) < 0 )){ // not <= because ties are covered above.
								///////////////////////////////////////////////////////////////////////////////////////////////
								// WE DO NOT HAVE VOTES YET!!!, so continue until we hit the slipScore basement.
								///////////////////////////////////////////////////////////////////////////////////////////////
								continue;

							} else {
								///////////////////////////////////////////////////////////////////////////////////////////////
								// WE DO NOT HAVE VOTES, BUT WE SLIPPED PAST THE CUTOFF SPECIFIED BY THE USER.
								///////////////////////////////////////////////////////////////////////////////////////////////
								TaxaLogger.logDebug("Analyzing OTU:" + otuId + " and gi:" + gi + " - Although we considered the top: " + votes.size() + " elements, we slipped past the slip-score.\n", null );
								break;
							}
						}
					}
				}

				String[][] taxonomicVotes = null; 
				boolean reportUnculturedCritters = false;

				if( votes.size() > 0 ){
					// Real votes.
					taxonomicVotes = new String[votes.size()][8];
					votes.toArray(taxonomicVotes);
					reportUnculturedCritters = false;

				} else if( votes.size() == 0 & candidates.size() > 0 ){
					// The votes at the top were uncultured; still, report a RESTRICTED taxonomic status.
					// We need to be able to 
					taxonomicVotes = new String[votesForARainyDay.size()][8];
					votesForARainyDay.toArray(taxonomicVotes);
					reportUnculturedCritters = true;
				}

				this.summarizeTaxonomicInfoForTargetSequence(writer, taxonomicVotes, taxonomicInfo, otuId, maximumBitScore, bitScore, eValue, reportUnculturedCritters);
			}

			for( String key : this.blastResultSequences.keySet()){
				if( !this.blastResults.containsKey( key )){
					// We have identified a sequence that has no significant similarity in the 
					// database under consideration (i.e. results were NOT written to file after mpiBLAST)
					this.summarizeTaxonomicInfoForTargetSequence(writer, null, null, key, -1, -1, -1, false);
				}
			}

			return true;

		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
			TaxaLogger.logError( "Error during vote: " + ioEx.getMessage(), null );

		} finally {
			if( writer != null ){
				writer.flush();
				writer.close();
			}
		}

		return false;
	}

	/*
	 * agreement thresholds for taxonomic status...
	 * THIS SECTION SHOULD BE REWRITTEN TO FORCEABLY ALIGN AGAINST RANK ENUM
	 */
	double[] thresholds =  new double[]{ 1, 1, 1, 1, 1, (double)2/3 };
	private double[] getThresholds(){
		return thresholds;
	}

	public void setThresholds( String cutoffs ){
		String[] cuts = (cutoffs.split(","));

		this.thresholds = new double[8];
		for( int i = 0; i < cuts.length; i++ ){
			if( cuts[i].indexOf("/") != -1 ){
				String[] dividing = cuts[i].split("/");
				double num = Double.parseDouble( dividing[0] );
				double denom = Double.parseDouble( dividing[1] );
				thresholds[i] = (double)num/denom;

			} else {
				thresholds[i] = Double.parseDouble(cuts[i]);
			}
		}
	}

	public void setThresholds( double[] cutoffs ){
		if( cutoffs.length == 6 ){
			this.thresholds = cutoffs;
		}
	}

	private void summarizeTaxonomicInfoForTargetSequence( BufferedWriter writer, String[][] ballots, Hashtable<String, Taxonomy> taxonomicInformation, 
			String otuId, double maximumBitScoreForTarget, double winningLevelBitScore, double winningEvalue, boolean uncultured ){

		// Exit early if there's nothing to report. :)
		if((null == ballots ) || ( null != ballots & ballots.length == 0 )){
			//omitted everything under consideration b/c culturables.size() == 0 at the 
			//scored-scale of interest (that is to say: if NOT USING a minscore, then, iterating
			//over all of the results yielded no cultured critters).
	        this.writeResult(writer, otuId, null, null, null, null, String.valueOf(winningEvalue), String.valueOf(winningLevelBitScore), 
	        		String.valueOf(maximumBitScoreForTarget), null, null, null, 0, null);
			TaxaLogger.logDebug( "There are no (valid) blast candidates for target-sequence: " + otuId, null );
	        return;
		}

		// 1.) Transpose the taxonomic ranks across the candidates...
		String[][] tmp = new String[8][ballots.length];
		for( int j = 0; j < ballots.length; j++ ){
			for( int k = 0; k < ballots[j].length; k++ ){
				tmp[k][j] = ballots[j][k];
			}

		} ballots = tmp; // transposed.

		// 2.) Iterate over taxonomic ranks; start at species and work UP to allow us to (hoffentlich) -succeed early-
		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		// VOTING-ALORITHM (votes b/w ties) at the spp, genus, family levels AND UP...
		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		int levelOfRank;
		for( Rank rank : Rank.values()){

			// Outputs in order given; here, we work up from SPP, GENUS, FAM, ..., PHYLUM...
			// This numeric code allows us to treat our Rank - enum as an int (for array management).
			levelOfRank = rank.getCode();

			// Determine the taxonomic assignments for this level (e.g. if taxonomic rank is SPECIES, this pulls out all of the votes at the SPECIES level).
			String[] votesAtThisRank = ballots[levelOfRank];
			Hashtable<String, Integer> counts = new Hashtable<String, Integer>();

			for( int vote_i = 0; vote_i < votesAtThisRank.length; vote_i++ ){ // length == number of voters;

				String currentVote = votesAtThisRank[vote_i];

				if( currentVote == null ){
					System.err.println("EXITING NOW. Null Entry From Taxonomy Database.");
					System.err.println("Using element: " + vote_i );
					System.err.println("i.e. taxonId (numeric): " + ballots[levelOfRank][vote_i]);
					System.err.println("key: " + currentVote);
					System.err.println("Rank-level: " + rank.toString());
					System.err.println("Folder: " + "" );
					System.exit(1);
				}

				// if the currently analyzed read didn't vote at this level (TODO: maybe the call is a SPUH) then ignore it (at this level of voting).
				if( "0".equals( currentVote )){ // we report 0s in Taxonomy.java for null Taxon-Ids.
					continue;
				}

				Integer count;
				if( counts.containsKey( currentVote )){
					// because we vote over taxon-ids (and assign literal-string names later), this works whether or not there exists a strain-id name concatenated to the spp name.
					count = counts.get( currentVote ); // e.g. Olpidium (gen) or Olpidium brassicae (spp) votes.

				} else {
					count = new Integer(0);
				}

				count++;
				counts.put( currentVote , count );
			} // At the conclusion of this loop, the votes will have been tallied. Next: check below to see if any taxonomic designation 'won'.

			String taxonId = null;
			Enumeration<String> e = counts.keys(); 

			while( e.hasMoreElements()){

				// Get the taxonomic-id and its count.
				taxonId = e.nextElement(); // retrieve the taxon-vote and the number of votes in its favor. 
				double votingThresholdTaxonomicRank_i = (double)getThresholds()[levelOfRank];
				double proportion = counts.get( taxonId ).doubleValue() / (double) votesAtThisRank.length;

				if( Double.compare( proportion, votingThresholdTaxonomicRank_i ) >= 0 ){ // does the vote count exceed our threshold (above)?
					break;
				}

				taxonId = null;
			}

			// We already have the taxonomic information (used it for voting) so now use it to describe the winner @ the winning rank.
			if( taxonId != null ){

				Taxonomy taxonomy_i = null;
				BlastCandidate candidate_i = null;
				ArrayList<BlastCandidate> candidates = this.blastResults.get( otuId );

				for( int i = 0; i <= candidates.size(); i++ ){ // iterate over the candidates, to list the GIs, and to determine if this GI has the winning taxonId.
					candidate_i = candidates.get(i);
					taxonomy_i = taxonomicInformation.get( candidate_i.getId()); // stored by GI, but we also need candidate information.
					if( taxonomy_i.getTaxonIds()[levelOfRank] == taxonId ){ // this asks, if, at this rank, the winning TaxonId is for this taxon.
						break;
					}
				}

				String taxonomicAffinity = null;
				switch( rank ){
				case GENUS:
					taxonomicAffinity = taxonomy_i.getGenus();     
					break;

				case FAMILY:
					taxonomicAffinity = taxonomy_i.getFamily();     
					break;

				case ORDER:
					taxonomicAffinity = taxonomy_i.getOrder();     
					break;

				case CLASS:
					taxonomicAffinity = taxonomy_i.getTaxonClass();     
					break;

				case PHYLUM:
					taxonomicAffinity = taxonomy_i.getPhylum();     
					break;

				case KINGDOM:
					taxonomicAffinity = taxonomy_i.getKingdom();     
					break;
				}

				// what if multiple GIs are associated with the same taxonomic-Id? that is why we 
				// take the winning values from the calling function, instead of from the matching candidate (whose GI is associated with the right taxonId)
				this.writeResult( 
						writer,
						otuId,
						candidate_i.getId(), // representative GI (there are bound to be multiple GIs at the same 'score', but this one is a suitable example).
						taxonId, 
						rank.toString(), 
						( uncultured ? taxonomicAffinity + " - Uncultured" : taxonomicAffinity ),
						( uncultured ? "NA" : String.valueOf( winningEvalue )),
						( uncultured ? "NA" : String.valueOf( winningLevelBitScore )), 
						String.valueOf( maximumBitScoreForTarget ), 
						candidate_i.getPercentageIdentity(), 
						format(((double) candidate_i.getLength()) / (Integer)this.blastResultSequenceLength.get( otuId ), true), 
						( uncultured ? taxonomy_i.getFullTaxonomy( rank ) + "Uncultured" : taxonomy_i.getFullTaxonomy( rank )), 
						counts.get(taxonId).intValue(), 
						new StringBuffer( String.valueOf( counts.get(taxonId).intValue())).append(":").append( votesAtThisRank.length ).toString());
				return;
			}
		}

        this.writeResult(writer, otuId, null, null, null, "No agreement at Kingdom", String.valueOf(winningEvalue), String.valueOf(winningLevelBitScore), 
        		String.valueOf(maximumBitScoreForTarget), null, null, null, 0, null);
		TaxaLogger.logDebug( "No agreement at the kingdom level for: " + otuId, null );
	}

	java.text.DecimalFormat df = new java.text.DecimalFormat("#.##");
	private String format( double percent, boolean truncateTo1 ){
		if( Double.compare((double)1, percent) < 0 ){
			return "1.0";

		} else {
			return( df.format( percent ));
		}
	}

	private boolean writeResult( BufferedWriter writer, String fastaId, String gi, String taxonid, String level, String winningTaxon, 
					String evalue, String winningScore, String topScore, String percentIdentity, 
					String coverage, String fullTaxonomy, int votes, String votesOutOf ){

		try{
	        writer.write( fastaId );
	        writer.write('\t');
	        writer.write(( gi != null ? gi.trim() : "NA"));
	        writer.write('\t');
	        writer.write(( taxonid != null ? taxonid.trim() : "NA"));
	        writer.write('\t');
	        writer.write(( level != null ? level.trim() : "NA"));
	        writer.write('\t');
	        writer.write(( winningTaxon != null ? winningTaxon.trim() : "Uncertain"));
	        writer.write('\t');
	        writer.write(( evalue != null ? evalue.trim() : "NA" ));
	        writer.write('\t');
	        writer.write(( winningScore != null ? winningScore.trim() : "NA" ));
	        writer.write('\t');
	        writer.write(( topScore != null ? topScore.trim() : "NA" ));
	        writer.write('\t');
	        writer.write(( percentIdentity != null ? percentIdentity.trim() : "NA" ));
	        writer.write('\t');
	        writer.write(( coverage != null ? coverage : "NA" ));
	        writer.write('\t');
	        writer.write(( fullTaxonomy != null ? fullTaxonomy.trim() : "NA" ));
	        writer.write('\t');
	        writer.write(( votes == 0 ? "NA" : String.valueOf( votes )));
	        writer.write('\t');
	        writer.write(( votesOutOf == null ? "NA" : String.valueOf( votesOutOf )));
	        writer.write('\n');
	        return true;

		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
			TaxaLogger.logError( "Error while announcing the winner: " + ioEx.getMessage(), fastaId );

		} finally {}

		return false;
	}

	//	tally.vote( blastResultsFileName, votingGroup, slippageTolerance, taxonomicCutoffs, percentIdentity );
	protected boolean vote( String blastResultsFileName, String votingGroup, double slippageTolerance, String cutoffs, double percentIdentity ){

		try {
			// 1.) build the hashmap of arrayListed blast-candidates.
			
			this.parseMPIblastOutput( blastResultsFileName );

			// 2.) filter down to 80% overlap
			this.filterResultsOnTargetLength( blastResultsFileName.replace( "_out", "" ));

			// 3.) find the maximum bit scores
			this.determineTopBitScores();

			// 4.) determine the winner.
			this.identifyTheWinners( blastResultsFileName.replaceAll("_out", ""), votingGroup, slippageTolerance, cutoffs, percentIdentity );

		} catch (IOException ioEx ){
			
		} finally {}

		return false;
	}
}

//private HashMap<String, String> originalSequences = null;
//private HashMap<String, Integer> originalSequenceLengths = null;