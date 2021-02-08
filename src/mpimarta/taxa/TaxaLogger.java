package mpimarta.taxa;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;

import org.apache.log4j.Logger;

public final class TaxaLogger {

	private static SimpleDateFormat fmt;
	private static Logger log4j = Logger.getLogger("mpimarta.taxa");

	private synchronized static String getTime(){
		fmt = new SimpleDateFormat();
		fmt.setCalendar(Calendar.getInstance());
		fmt.applyPattern("MMddyy-HHmmss");
		return fmt.format( fmt.getCalendar().getTime());
	}

	public synchronized static void logDebug( String debug, String shortid ){
		log4j.debug( new StringBuffer(( shortid == null ? "" : "Id: ".concat( shortid ))).append( " time: " ).append( 
				getTime()).append( " " ).append( debug ).toString());
	}

	public synchronized static void logError( String error, String shortid ){
		log4j.error( new StringBuffer(( shortid == null ? "" : "Id: ".concat( shortid ))).append( " time: " ).append( 
				getTime()).append( " " ).append( error ).toString());
	}

	public synchronized static void logInfo( String info, String shortid ){
		log4j.info( new StringBuffer(( shortid == null ? "" : "Id: ".concat( shortid ))).append( " time: " ).append( 
				getTime()).append( " " ).append( info ).toString());
	}

	public static void updateVotingHistory( String outputdir, String groupId, String options ){

		BufferedWriter writer = null;

		try{

	        writer = new BufferedWriter(
	        		new OutputStreamWriter( new FileOutputStream( 
	        				new StringBuffer( outputdir ).append( "/history.txt" ).toString(), true ), "UTF-8"));

	        writer.write("----------------------------------------------------------------\n");
	        writer.write("groupId: ");
	        writer.write(groupId);
	        writer.write('\n');
	        writer.write( options );
	        writer.write('\n');
	        writer.flush();
	        writer.close();

		} catch( IOException ioEx ){
			ioEx.printStackTrace(System.err);
			log4j.error( "Error while recording the history.");
			log4j.error( ioEx.getMessage());

		} finally {}
	}
}
