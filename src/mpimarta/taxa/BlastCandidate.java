package mpimarta.taxa;

public class BlastCandidate implements Comparable<BlastCandidate> {

	private final int BLAST_COLUMN_ID = 2;
	private final int BLAST_COLUMN_PERC_IDENTITY = 3;
	private final int BLAST_COLUMN_ALIGN_LENGTH = 4;
	private final int BLAST_COLUMN_MINIMUM_EVALUE = 11;
	private final int BLAST_COLUMN_BITSCORE = 12;

	private String id;
	private String fastaId;
	private String percentageIdentity;
	private int length;
	private double minimumEvalue;
	private double bitscore;

	public BlastCandidate( String fastaId, String[] fields ){

		this.fastaId = fastaId;

		String[] ids = fields[BLAST_COLUMN_ID-1].split("\\|")[1].split("\\."); //[0]; // pipes separate the ids, but there are often dots inside the names THAT ARE NOT INCLUDED IN THE ALIGNMENT NAMES (to denote strand). e.g. emb|AJ242887.1|RSP242887
		this.id = ids[0];
		this.percentageIdentity = fields[BLAST_COLUMN_PERC_IDENTITY - 1];
		this.length = Integer.parseInt( fields[BLAST_COLUMN_ALIGN_LENGTH - 1]);

// minimum E-value and bitscore.
		String tmp = fields[BLAST_COLUMN_MINIMUM_EVALUE - 1];
		if( tmp.startsWith("e")){ tmp = "1".concat(tmp); }
		this.minimumEvalue = Double.parseDouble(tmp );
		this.bitscore = Double.parseDouble( fields[ BLAST_COLUMN_BITSCORE - 1] );
	}

//	gi
	public String getId() {
		return this.id;
	}

	public void setId(String id) {
		this.id = id;
	}

//	short-id
	public String getFastaId(){
		return this.fastaId;
	}

	public void setFastaId( String fastaId ){
		this.fastaId = fastaId;
	}

//	perc-identity (%age similarity b/w target and this hit for the given overlap);
	public String getPercentageIdentity(){
		return this.percentageIdentity;
	}

//	length itself.
	public int getLength(){
		return this.length;
	}

//	minimum e-value
	public double getMinimumEvalue(){
		return this.minimumEvalue;
	}

//	bitscore
	public double getBitscore(){
		return this.bitscore;
	}

	@Override 
	public int compareTo( BlastCandidate o ){
		// flip the order of the two objects to reverse sort.
		return Double.compare(o.getBitscore(), this.getBitscore());
	}
}
