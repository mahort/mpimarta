package mpimarta.taxa;

import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

public enum Rank {

	// load backwards for the iterator (which works 'upwards' from spp to phylum --so far--).
	GENUS(5), FAMILY(4), ORDER(3), CLASS(2), PHYLUM(1), KINGDOM(0);

	private static final Map<Integer,Rank> lookup = new HashMap<Integer,Rank>();

	static {
		for(Rank s : EnumSet.allOf(Rank.class)){
			lookup.put(s.getCode(), s);
		}
	}

	private int code;

	private Rank(int code){
	     this.code = code;
	}

	public int getCode() { return code; }

	public static Rank get( int code ){ 
		return lookup.get(code); 
	}
}