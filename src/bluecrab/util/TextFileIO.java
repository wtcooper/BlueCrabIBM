package bluecrab.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;

/**
 * Utility file input/output class with synchronized methods to print or read lines.  
 *  
 * @author Wade.Cooper
 *
 */

public class TextFileIO {

	private String filename; 
	private PrintWriter out= null; 
	private BufferedReader reader = null; 
	private boolean appendFile = false;
	private String delim="\t"; //tab-delimited is default

	public TextFileIO(String filename){
		this.filename = new File(filename).getAbsolutePath();  
	}




	/**
	 * Prints a string via synchronized method. If the argument is null then the string "null" 
	 * is printed.  Otherwise, the string's characters are converted into bytes according to the 
	 * platform's default character encoding, and these bytes are written in exactly 
	 * the manner of the write(int) method.
	 * 
	 * @param s The String to be printed
	 */
	public synchronized void print(String s){
		if (out == null) out = getWriter();
		out.print(s);
	}

	/**
	 * Prints a String via syncrhonized method and then terminates the line. 
	 * This method behaves as though it invokes print(String) and then println().
	 * 
	 * @param s
	 */
	public synchronized void println(String s){
		if (out == null) out = getWriter();
		out.println(s);
	}

	/**
	 * Prints all Strings or numeric values, separated by commas, via syncrhonized method
	 *  and then terminates the line. 
	 * 
	 * 
	 * @param s
	 */
	public synchronized void println(Object... sArr){
		if (out == null) out = getWriter();
		for (int i=0; i<sArr.length; i++) {
			String s=null;
			//print a regular string
			if (sArr[i] instanceof String) {
				s=(String) sArr[i];
				//don't add tab to end of line, so need to check if last element or not
				if (i<sArr.length-1) out.print(s+delim);
				else out.println(s);

			}

			//if the object is a double array
			else if (sArr[i] instanceof double[] ) { //if a primitive array
				double[] temp = (double[]) sArr[i]; 
				//iterate over the array and put in tab separator
				for (int j=0; j<temp.length; j++){
					if (i<sArr.length-1) out.print(temp[j]+delim);
					else {
						if (j<temp.length-1) out.print(temp[j]+delim);
						else out.println(temp[j]);
					}
				}
			}

			//if the object is a double array
			else if (sArr[i] instanceof int[] ) { //if a primitive array
				int[] temp = (int[]) sArr[i]; 
				//iterate over the array and put in tab separator
				for (int j=0; j<temp.length; j++){
					if (i<sArr.length-1) out.print(temp[j]+delim);
					else {
						if (j<temp.length-1) out.print(temp[j]+delim);
						else out.println(temp[j]);
					}
				}
			}

			//if the object is a double array
			else if (sArr[i] instanceof long[] ) { //if a primitive array
				long[] temp = (long[]) sArr[i]; 
				//iterate over the array and put in tab separator
				for (int j=0; j<temp.length; j++){
					if (i<sArr.length-1) out.print(temp[j]+delim);
					else {
						if (j<temp.length-1) out.print(temp[j]+delim);
						else out.println(temp[j]);
					}
				}
			}

			//if the object is a double array
			else if (sArr[i] instanceof float[] ) { //if a primitive array
				float[] temp = (float[]) sArr[i]; 
				//iterate over the array and put in tab separator
				for (int j=0; j<temp.length; j++){
					if (i<sArr.length-1) out.print(temp[j]+delim);
					else {
						if (j<temp.length-1) out.print(temp[j]+delim);
						else out.println(temp[j]);
					}
				}
			}

			
			//capture anything else and convert to string representation
			else {
				s=String.valueOf(sArr[i]);
				//don't add tab to end of line, so need to check if last element or not
				if (i<sArr.length-1) out.print(s+delim);
				else out.println(s);
			}
		}
	}

	/**
	 * Prints a String via syncrhonized method and then terminates the line. 
	 * This method behaves as though it invokes print(String) and then println().
	 * 
	 * @param s
	 */
	public synchronized void println(){
		if (out == null) out = getWriter();
		out.println();
	}


	/**
	 * Reads a line of text. A line is considered to be terminated by any one of a line feed ('\n'), 
	 * a carriage return ('\r'), or a carriage return followed immediately by a linefeed.
	 * 
	 * @return A String containing the contents of the line, not including any line-termination 
	 * characters, or null if the end of the stream has been reached
	 */
	public synchronized String readLine(){
		if (reader == null) reader = getReader();
		try {
			return reader.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}




	/**
	 * Close up shop when done.
	 */
	public void close() {
		if (out != null) out.close();
		if (reader != null) {
			try { 
				reader.close();
			} catch (IOException e) {e.printStackTrace();
			}
		}
	}


	/**
	 * Sets whether the writer should append the file or replace the file.  Default is to replace.
	 * 
	 * @param append
	 */
	public void setAppend(boolean append) {
		this.appendFile = append;
	}

	/**
	 * Get a PrintWriter.  
	 * 
	 * @return
	 */
	private PrintWriter getWriter(){

		//erase any old files if writing to file and not appending
		if (!appendFile) new File(filename).delete();
		File fFile = new File(filename);

		try { 

			out= new PrintWriter(new FileWriter(fFile, true));
			return out;
		} catch (IOException e) {e.printStackTrace();
		}
		return null;
	}



	/**
	 * Get a BufferedReader.  Use as:
	 * BufferedReader reader = textFileIO.getReader();
	 * 
	 * 			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
	 *				String[] tokens = line.split("\t");
	 *				...
	 *			} 
	 * 
	 * @return
	 */
	private BufferedReader getReader() {

		File fFile = new File(filename);
		try { 
			reader = new BufferedReader(new FileReader(fFile));
			return reader; 
		} catch (IOException e) {e.printStackTrace();
		}
		return null;
	}



	/**
	 * Copy a file from source to destination.
	 * 
	 * @param srFile
	 * @param dtFile
	 */
	public static void copyfile(String srFile, String dtFile){
		try{
			File f1 = new File(srFile);
			File f2 = new File(dtFile);
			InputStream in = new FileInputStream(f1);

			//For Append the file.
			//  OutputStream out = new FileOutputStream(f2,true);

			//For Overwrite the file.
			OutputStream out = new FileOutputStream(f2);

			byte[] buf = new byte[1024];
			int len;
			while ((len = in.read(buf)) > 0){
				out.write(buf, 0, len);
			}
			in.close();
			out.close();
			System.out.println("File copied.");
		}
		catch(FileNotFoundException ex){
			System.out.println(ex.getMessage() + " in the specified directory.");
			System.exit(0);
		}
		catch(IOException e){
			System.out.println(e.getMessage());  
		}
	}



	/**
	 * Set the delimiter in the file (default is tab-delimited).
	 * 
	 * @param delim
	 */
	public void setDelim(String delim) {
		this.delim = delim;
	}



}
