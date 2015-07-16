package bluecrab;

import java.util.ArrayList;
import java.util.HashMap;

import bluecrab.util.TextFileIO;



public final class BlueCrabData {

	


	/**
	 * Returns the FL temperature data
	 * @return
	 */
	public static double[][] getTempData() {
		
		double[][] tempData = new double[538][2]; //column1 = age, column2=temp (deg C) 

		TextFileIO tempDataFile = new TextFileIO("input/temp.txt");
		int counter = 0;
		for (String line = tempDataFile.readLine(); line != null; line = tempDataFile.readLine()) { 
			String[] tokens = line.split("\t"); 
			tempData[counter][0] = Double.parseDouble(tokens[0]);
			tempData[counter][1] = Double.parseDouble(tokens[1]);
			counter++;
		}

		tempDataFile.close();

		return tempData;
	}
	
	
	/**
	 * Returns the FL average empirical temperature data
	 * @return
	 */
	public static double[] getSimTempData(String file) {
		
		double[] tempData = new double[366]; //column1 = age, column2=temp (deg C) 

		TextFileIO tempDataFile;  
		tempDataFile= new TextFileIO(file);
		
		tempDataFile.readLine(); //read in header
		
		int counter = 0;
		for (String line = tempDataFile.readLine(); line != null; line = tempDataFile.readLine()) { 
			String[] tokens = line.split("\t"); 
			tempData[counter] = Double.parseDouble(tokens[1]);
			counter++;
		}

		tempDataFile.close();

		return tempData;
	}

	
	/**
	 * Returns the FL temperature data
	 * @return
	 */
	public static double[] getSimOvigData(String file) {
		
		double[] tempData = new double[12]; //column1 = age, column2=temp (deg C) 

		TextFileIO tempDataFile;  
		tempDataFile= new TextFileIO(file);
		
		tempDataFile.readLine(); //read in header
		
		int counter = 0;
		for (String line = tempDataFile.readLine(); line != null; line = tempDataFile.readLine()) { 
			String[] tokens = line.split("\t"); 
			tempData[counter] = Double.parseDouble(tokens[1]);
			counter++;
		}

		tempDataFile.close();

		return tempData;
	}
	
	/**
	 * Return the Mississippi pond growth temperature data.
	 * 
	 * @return
	 */
	public static ArrayList<double[]> getMSGrowthData() {
		ArrayList<double[]> data = new ArrayList<double[]>();
		TextFileIO tempDataFile = new TextFileIO("input/MSData.txt");
		//NumObs	133	258	121	14	169	247	127
		double[] data1 = new double[133];
		double[] data2 = new double[258];
		double[] data3 = new double[121];
		double[] data4 = new double[14];
		double[] data5 = new double[169];
		double[] data6 = new double[247];
		double[] data7 = new double[127];

		tempDataFile.readLine(); //read the header
		int counter=0;
		for (String line = tempDataFile.readLine(); line != null; line = tempDataFile.readLine()) { 
			String[] tokens = line.split("\t"); 
			double temp1 = Double.parseDouble(tokens[1]);
			double temp2 = Double.parseDouble(tokens[2]);
			double temp3 = Double.parseDouble(tokens[3]);
			double temp4 = Double.parseDouble(tokens[4]);
			double temp5 = Double.parseDouble(tokens[5]);
			double temp6 = Double.parseDouble(tokens[6]);
			double temp7 = Double.parseDouble(tokens[7]);
			
			if (temp1 != 999) data1[counter] = temp1;
			if (temp2 != 999) data2[counter] = temp2;
			if (temp3 != 999) data3[counter] = temp3;
			if (temp4 != 999) data4[counter] = temp4;
			if (temp5 != 999) data5[counter] = temp5;
			if (temp6 != 999) data6[counter] = temp6;
			if (temp7 != 999) data7[counter] = temp7;
			
			counter++;
		}
		data.add(data1);
		data.add(data2);
		data.add(data3);
		data.add(data4);
		data.add(data5);
		data.add(data6);
		data.add(data7);
		
		return data;
	}
	
	
	public static ArrayList<double[][]> getAvgGrowthData() {
		
		ArrayList<double[][]> data = new ArrayList<double[][]>();

		double[][] fData = new double[58][4];
		double[][] mData = new double[58][4];

		//########################
		//read temp data
		TextFileIO growthFile = new TextFileIO("input/avgGrowth.txt");

		int counter = 0;
		
		for (String line = growthFile.readLine(); line != null; line = growthFile.readLine()) { 
			if (line.isEmpty() || line.contains("age")) continue; //ignore lines with comments

			String[] tokens = line.split("\t"); 

			fData[counter][0] = Double.parseDouble(tokens[0]);
			fData[counter][1] = Double.parseDouble(tokens[1]);
			fData[counter][2] = Double.parseDouble(tokens[2]);
			fData[counter][3] = Double.parseDouble(tokens[3]);

			mData[counter][0] = Double.parseDouble(tokens[0]);
			mData[counter][1] = Double.parseDouble(tokens[4]);
			mData[counter][2] = Double.parseDouble(tokens[5]);
			mData[counter][3] = Double.parseDouble(tokens[6]);

			counter++;

		}
		growthFile.close();
		
		data.add(mData);
		data.add(fData);
		return data;

	}
		

	/**
	 * Return the male and female specific data, males first and females second in ArrayList
	 * @return
	 */
	public static ArrayList<double[][]> getIndGrowthData() {
		ArrayList<double[][]> data = new ArrayList<double[][]>();

		//### FEMALE Data ###
		int counter = 0;
		double[][] fData = new double[857][2];
		TextFileIO fGrowthFile = new TextFileIO("input/FemalePondGrowth.txt");
		for (String line = fGrowthFile.readLine(); line != null; line = fGrowthFile.readLine()) { 
			if (line.isEmpty() || line.contains("age")) continue; //ignore lines with comments
			String[] tokens = line.split("\t"); 
			fData[counter][0] = Double.parseDouble(tokens[1]);		//age
			fData[counter][1] = Double.parseDouble(tokens[3]);		//size
			counter++;

		}
		fGrowthFile.close();


		//### MALE Data ###
		counter = 0;
		double[][] mData = new double[1049][2];
		TextFileIO mGrowthFile = new TextFileIO("input/MalePondGrowth.txt");
		for (String line = mGrowthFile.readLine(); line != null; line = mGrowthFile.readLine()) { 
			if (line.isEmpty() || line.contains("age")) continue; //ignore lines with comments
			String[] tokens = line.split("\t"); 
			mData[counter][0] = Double.parseDouble(tokens[1]);		//age
			mData[counter][1] = Double.parseDouble(tokens[3]);		//size
			counter++;
		}
		mGrowthFile.close();

		
		data.add(mData);
		data.add(fData);
		return data;
	}

	
	/**
	 * Return the male and female specific data for MS as a HashMap, where map key
	 * is a String of Sex+PondNum, e.g., M1 for males in the first pond.  ArrayList of each 
	 * sex and pond is just the CW (no need for age since all same age)
	 * 
	 * @return
	 */
	public static HashMap<String,ArrayList<Double>> getMSIndGrowthData() {
		HashMap<String,ArrayList<Double>> data = new HashMap<String,ArrayList<Double>>();

		TextFileIO growthFile = new TextFileIO("input/MSSizeData.txt");
		for (String line = growthFile.readLine(); line != null; line = growthFile.readLine()) { 
			if (line.isEmpty() || line.contains("Trial")) continue; //ignore header file
			String[] tokens = line.split("\t"); 
			String key;
			if (tokens[1].equals("1")) key="M"+tokens[0];
			else key="F"+tokens[0];
			ArrayList<Double> list = data.get(key);
			if (list==null) {
				list=new ArrayList<Double>();
				data.put(key, list);
			}
			list.add(Double.parseDouble(tokens[2]));
		}
		growthFile.close();

		return data;
	}
	
	
	 /**
	  * Gets a double array [obs][0=age,1=size] for the full growth dataset
	  * @return
	  */
	public static double[][] getIndGrowthDataBothSexes() {

		//### FEMALE Data ###
		int counter = 0;
		double[][] data = new double[1714][2];
		TextFileIO growthFile = new TextFileIO("input/PondGrowth.txt");
		for (String line = growthFile.readLine(); line != null; line = growthFile.readLine()) { 
			if (line.isEmpty() || line.contains("age")) continue; //ignore lines with comments
			String[] tokens = line.split("\t"); 
			data[counter][0] = Double.parseDouble(tokens[1]);		//age
			data[counter][1] = Double.parseDouble(tokens[3]);		//size
			counter++;

		}
		growthFile.close();
		
		return data;
	}
	
	
	 /**
	  * Gets a double array [obs][0=age,1=size] for the full growth dataset
	  * @return
	  */
	public static double[][][] getEstMeanSizes(String sexIn) {

		int maleCounter= 0;
		int femaleCounter=0;
		double[][][] data = new double[538][2][2];	//[day observations][sex][day,size]
		TextFileIO sizesFile = new TextFileIO("output/blueCrabAvgSizes_"+sexIn+".txt");

		for (String line = sizesFile.readLine(); line != null; line = sizesFile.readLine()) { 
			if (line.isEmpty() || line.contains("age")) continue; //ignore lines with comments
			String[] tokens = line.split("\t"); 
			int sex=Integer.parseInt(tokens[0]);
			if (sex==0) {
				data[maleCounter][0][0] = Double.parseDouble(tokens[1]);		//age
				data[maleCounter][0][1] = Double.parseDouble(tokens[2]);		//size
				maleCounter++;
			}
			else {
				data[femaleCounter][1][0] = Double.parseDouble(tokens[1]);		//age
				data[femaleCounter][1][1] = Double.parseDouble(tokens[2]);		//size
				femaleCounter++;
			}

		}
		sizesFile.close();
		
		return data;
	}



	
	
	/**
	 * Returns mean GPM for the sizes 1-200mm.
	 * 
	 * @return
	 */
	public static double[] getEstMeanGPM(String sexIn) {
		double[] data = new double[200];
		TextFileIO gpmFile = new TextFileIO("output/meanGPM_" + sexIn + ".txt");
		for (String line = gpmFile.readLine(); line != null; line = gpmFile.readLine()) { 
			if (line.contains("CW")) continue; //ignores header
			String[] tokens = line.split("\t"); 
			int cw=Integer.parseInt(tokens[0]);
			double gpm = Double.parseDouble(tokens[1]);
			data[cw]=gpm;
		}
		
		return data;
	}
	
}
