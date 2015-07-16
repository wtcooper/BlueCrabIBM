package bluecrab.ecj;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.StatUtils;

import sim.util.distribution.Beta;
import sim.util.distribution.Normal;
import sim.util.distribution.Uniform;
import bluecrab.BlueCrab;
import bluecrab.BlueCrabData;
import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.simple.SimpleFitness;
import ec.simple.SimpleProblemForm;
import ec.util.MersenneTwisterFast;
import ec.util.Parameter;
import ec.vector.DoubleVectorIndividual;

/**
 * This fits the full suite of paramters (including variability params) to both Florida and
 * MS pond growth data, using a molt process temperature dependent IBM. 
 * 
 * Here, model crabs for each of 8 ponds (1 FL, 7 MS) using the measured temperature time series
 * data for each pond.  For each time step when crabs were measured (~ every week for FL,
 * end time for MS), compute the ave+/-SD for the modeled crabs.  Use this expected 
 * mu and sigma in the full log likelihood to compare to the individual observed sizes.  The sigma will
 * fit the variance of the length at age appropriately, thus allowing the appropriate selection of 
 * GPM_stDev and the shifted exponential alpha parameters.  
 * 
 * @author Wade.Cooper
 *
 */


public class BlueCrabOptimize extends Problem implements SimpleProblemForm {

	private static final double idealFitness=1.0; 
	
	private static final boolean runOnlyFL = false;  //runs only the Florida data

	private static final long serialVersionUID = 1L;
	private static final Object LOCK_1 = new Object() {};
	public static final String P_NUMCRABS = "numCrabs";
	public static final double L50_MATURE = 125.202671;
	public static final double STD_MATURE = 2.909795;
	public static final double SEX_RATIO = .5;			//males:females
	public static final double SIZE_AT_SETTLEMENT = 2.5;
	public static final double STD_SIZE_AT_SETTLEMENT = .1;

	public static final double lowBounds=1.025;
	public static final double upBounds=1.5;
	public static final double maxXVal=200;

	public static int gpmErrType=0;	//0=lognormal, 1=beta, 2=normal

	public static int numCrabs; 	//read in from param file		

	public static int numTimeSteps = 538;		//For FL, start on 7/23/2010, day of settlement (62d prior to growth data)



	//######### FLORIDA ###########
	public static double[][] tempData;		//temperature data
	public static double[][] mObs;				//male observed data
	public static double[][] fObs;				//female observed data
	public static int[] fObsDays;					//age in days when observed 
	public static int[] mObsDays;					//age in days when observed 


	//######### MS ###########
	//Start days: 5/12/2009	8/26/2009	1/28/2010	6/2/2010	11/19/2010	5/26/2011	7/14/2011
	//Start and stopping CW's
	public static double[][] MSIniSize = {{23.71190476, 20.34619048, 15.53666667, 17.29416667, 15.081875, 18.55276382, 14.365},
		{8.1731736, 6.54891859, 3.404991745, 5.059549914, 3.247153991, 5.267741275, 5.073507405}};
	//number of time steps to step for MS data
	public static int[] MSTimeSteps = {133, 258, 121, 14, 169, 247, 127};
	public static HashMap<String,ArrayList<Double>> MSObsSize;

	public static ArrayList<double[]> MSTempData ;

	//holds the recycled crabs
	public static ArrayList<ArrayList<BlueCrab>>recycledCrabs;


	/**
	 * Setup method for ECJ.
	 * Runs a single time for the optimization, so here's where to read in static data and 
	 * parameters from the .params file
	 */
	public void setup(final EvolutionState state, final Parameter base) { 
		super.setup(state,base);

		if (!runOnlyFL){
			//get the temperature data
			MSTempData = BlueCrabData.getMSGrowthData();
			//get the observed data
			MSObsSize=BlueCrabData.getMSIndGrowthData();
		}

		tempData = BlueCrabData.getTempData();
		ArrayList<double[][]> data = BlueCrabData.getIndGrowthData();
		mObs = data.get(0);
		fObs = data.get(1);

		ArrayList<Integer> femTemp = new ArrayList<Integer>(); 
		ArrayList<Integer> maleTemp = new ArrayList<Integer>(); 
		Integer[] temp=new Integer[1];
		//set the observed days that data available 
		for (int i=0; i<mObs.length; i++) //get the age data for binary sort later
			if (!maleTemp.contains((int) mObs[i][0])) maleTemp.add((int) mObs[i][0]);  
		for (int i=0; i<fObs.length; i++) 
			if (!femTemp.contains((int) fObs[i][0])) femTemp.add((int) fObs[i][0]);  

		//convert from Integer array to int array for memory
		mObsDays = ArrayUtils.toPrimitive(maleTemp.toArray(temp));
		fObsDays = ArrayUtils.toPrimitive(femTemp.toArray(temp));

		//read in the number of crabs from parameter file
		numCrabs = state.parameters.getInt(base.push(P_NUMCRABS),null,1);

	}


	/**
	 * Main stepping method ECJ
	 */
	@Override
	public void evaluate(EvolutionState state, Individual ind,
			int subpopulation, int threadnum) {


		//######################################
		//Non-global data for individual parameter runs
		//######################################

		MersenneTwisterFast m = state.random[threadnum]; 
		Uniform uniform = new Uniform(m); 
		Normal normal = new Normal(0,1,m); 
		Beta betaDist = new Beta(1,1,m);


		double[][] mEst = new double[58][(int) (numCrabs/2.)];		//estimated male size at age
		double[][] fEst = new double[58][(int) (numCrabs/2.)];		//estimated female size at age

		//holds the ending CW estimates for the MS growth data
		//first dim is the trial/pond num, next are the individual sizes
		double[][] mMSEst = new double[7][(int) (numCrabs/2.)];
		double[][] fMSEst = new double[7][(int) (numCrabs/2.)];

		double gpm=1, meanGPM, logMu, logStDev, cv, var, alpha, beta, adjMean, fitness = 0;

		DoubleVectorIndividual ind2 = (DoubleVectorIndividual) ind; //genetically-unique individual




		//######################################
		//Parameters
		//######################################

		//GPM Parameters (spline and st dev)
		double[] knotsX={ind2.genome[0],ind2.genome[2],ind2.genome[4],ind2.genome[6],ind2.genome[8]};
		double[] knotsY={ind2.genome[1],ind2.genome[3],ind2.genome[5],ind2.genome[7],ind2.genome[9]};
		double GPM_stDev = ind2.genome[10];	//standard deviation for GPM

		//IP parameters; gamma and beta
		double gamma_a = ind2.genome[11];
		double gamma_b= ind2.genome[12];
		double gamma_hi = 100; //ind2.genome[13];	//Bunnell and Miller default 2.11, set high to void

		double beta_a = ind2.genome[14];
		double beta_b = ind2.genome[15];


		//######################################
		//Constraints
		//######################################


		//Constraint for spline knots in order
		for (int i=0; i<knotsX.length-1; i++){
			if (knotsX[i]>knotsX[i+1]	) {
				((SimpleFitness)ind2.fitness).setFitness(state,(float) 0.0,fitness == idealFitness);
				ind2.evaluated = true;
				return;
			}
		}

		//Gamma a must be < beta a 
		if (gamma_a>beta_a) {
			((SimpleFitness)ind2.fitness).setFitness(state,(float) 0.0,fitness == idealFitness);
			ind2.evaluated = true;
			return;
		}

		//Set up the spline interpolation
		SplineInterpolator cs = new SplineInterpolator();
		PolynomialSplineFunction spline=cs.interpolate(knotsX, knotsY);
		double lastgpm=0; 

		//Constraint to make sure properly bounded
		for (int i=1; i<200; i++){
			gpm=1;
				//need to set upper bound for function, else will throw error
				if (i>maxXVal) gpm=spline.value(maxXVal);
				else gpm=spline.value(i);

			//check to make sure in bounds
			if (gpm<lowBounds || gpm >upBounds) {
				((SimpleFitness)ind2.fitness).setFitness(state,(float) 0.0,fitness == idealFitness);
				ind2.evaluated = true;
				return;
			}

			//check to make sure gpm doesn't increase at end where limited data
			if (i>100 && gpm>lastgpm){
				((SimpleFitness)ind2.fitness).setFitness(state,(float) 0.0,fitness == idealFitness);
				ind2.evaluated = true;
				return;
			}

			lastgpm=gpm;
		}





		//######################################
		//Instantiate the populations
		//######################################


		ArrayList<BlueCrab> crabs = getCrabs();

		for (int i=0; i<8; i++) {	//loop through each pond trial; 0=florida, 1->7=MS

			if (runOnlyFL && i>0) continue; 

			for (int j=0; j<numCrabs; j++){

				BlueCrab crab;
				if (crabs.size() <= (i*numCrabs+j)){		//linear indexing of crabs
					crab = new BlueCrab();
					crabs.add(crab);
				}
				else crab = crabs.get(i*numCrabs+j);		

				crab.id=i;	//mark each crab id to belong to the appropriate trial

				//initialize the FL crabs
				if (i==0) {
					crab.initialize();
					crab.sex = 0; 
					if (j >= numCrabs/2.) crab.sex = 1; //set 1/2 crabs male and other half female

					//make lognormal since can't drop < 0; see Lognormal.java for this derivation 
					logStDev=Math.sqrt(Math.log(1+(STD_SIZE_AT_SETTLEMENT*STD_SIZE_AT_SETTLEMENT)/(SIZE_AT_SETTLEMENT*SIZE_AT_SETTLEMENT)));
					logMu = Math.log(SIZE_AT_SETTLEMENT) - 0.5*(logStDev*logStDev); 
					crab.CW = Math.exp(normal.nextDouble(logMu,logStDev));		

					//Pull IP from shifted exponential
					crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
							*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
							+ gamma_a*Math.pow(gamma_b,crab.CW);

					while (crab.IP > gamma_hi *  beta_a*Math.pow(beta_b,crab.CW)	) {
						crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
								*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
								+ gamma_a*Math.pow(gamma_b,crab.CW);
					}

					//this is only used for females -- logistic random number
					crab.sizeMature = L50_MATURE - Math.log((1/uniform.nextDoubleFromTo(0, 1))-1)*STD_MATURE;
				} // end if FL check

				//initialize the MS crabs
				else {

					crab.ageInDays=92;
					crab.degreeDays=0;
					crab.numMolts=0;
					crab.sex=0;
					if (j >= numCrabs/2.) crab.sex = 1; //set 1/2 crabs male and other half female

					logStDev=Math.sqrt(Math.log(1+(MSIniSize[1][i-1]*MSIniSize[1][i-1])/(MSIniSize[0][i-1]*MSIniSize[0][i-1])));
					logMu = Math.log(MSIniSize[0][i-1]) - 0.5*(logStDev*logStDev); 
					crab.CW = Math.exp(normal.nextDouble(logMu,logStDev));		

					crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
							*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
							+ gamma_a*Math.pow(gamma_b,crab.CW);

					while (crab.IP > gamma_hi *  beta_a*Math.pow(beta_b,crab.CW)	) {
						crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
								*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
								+ gamma_a*Math.pow(gamma_b,crab.CW);
					}

					//this is only used for females
					crab.sizeMature = L50_MATURE-Math.log((1/uniform.nextDoubleFromTo(0, 1))-1)*STD_MATURE;

				} // MS check

			} // loop over num crabs
		} //loop over trials/ponds




		//##################################################
		//run the time steps for JUST Ryan's FL growth 
		// I could combine the code in the i-loop over pond trials, but i was 
		// lazy in prettying it up 
		//##################################################


		for (int j=0; j<numCrabs; j++) {

			BlueCrab crab = crabs.get(j);

			for (int k=0; k<numTimeSteps; k++) {


				//if male, or if female and immature, then check if it should grow
				if (crab.sex == 0 || (crab.sex == 1  && crab.CW < crab.sizeMature)) {

					//(1) accummulate degree days
					double degDays = tempData[k][1] - 8.9;

					if (degDays > 0) 
						crab.degreeDays += degDays;

					//(2) check if molt, degree days > IP
					if (crab.degreeDays > crab.IP) {


						//get the mean GPM for a crab of this size
						meanGPM=1;
							//need to set upper bound for function, else will throw error
							if (crab.CW>maxXVal) meanGPM=spline.value(maxXVal);
							else meanGPM=spline.value(crab.CW);
						//some rounding errors in func.value; need to make sure never drops below
						if (meanGPM<lowBounds) meanGPM=lowBounds+0.001;

						//Get a variable GPM around the mean spline, using either lognorm or beta:

						//Lognormal distributed GPM
						if (gpmErrType==0){
							logStDev=Math.sqrt(Math.log(1+(GPM_stDev*GPM_stDev)/((meanGPM-lowBounds)*(meanGPM-lowBounds))));
							logMu = Math.log(meanGPM-lowBounds) - 0.5*(logStDev*logStDev); 
							gpm = lowBounds+Math.exp(normal.nextDouble(logMu,logStDev));
						}

						//Beta distributed GPM
						else if (gpmErrType==1) {
							//Beta distribution calculation of alpha and beta from sample mean & var
							adjMean=(meanGPM-lowBounds)/(upBounds-lowBounds);
							cv=GPM_stDev/meanGPM;	//set up the CV
							var=(cv*adjMean)*(cv*adjMean);	//calculate variance
							//calculate the alpha and beta from the adjusted mu/var
							alpha =((1-adjMean)/var-1/adjMean)*(adjMean*adjMean);
							beta = alpha*(1/adjMean-1);								
							gpm=(upBounds-lowBounds)*betaDist.nextDouble(alpha, beta)+lowBounds;
						}
						else if (gpmErrType==2){
							gpm=normal.nextDouble(meanGPM, GPM_stDev);
						}
						
						crab.CW *= gpm;

						crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
								*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
								+ gamma_a*Math.pow(gamma_b,crab.CW);

						while (crab.IP > gamma_hi *  beta_a*Math.pow(beta_b,crab.CW)	) {
							crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
									*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
									+ gamma_a*Math.pow(gamma_b,crab.CW);
						}

						crab.degreeDays = 0;  // reset to 0
					} // end check if crabs molt

				}//check if male or immature female in order to grow

				crab.ageInDays++;

				/*get the array index of the day measurement in observed data,
				 * and if there were observed data taken on this day, then store
				 * the predicted value for calculating fitness later on 
				 */
				int index=0;
				if (crab.sex == 0) index = Arrays.binarySearch(mObsDays, k+30);
				else index = Arrays.binarySearch(fObsDays, k+30);
				if ( index >= 0) {	//if positive, then is a match

					if (crab.sex == 0)
						mEst[index][j] = crab.CW;
					else 
						fEst[index][j-(int) (numCrabs/2.)] = crab.CW;	//females are in 2nd half of array
				}

			}// loop over time steps
		}// loop over crabs


		//##################################################
		//2nd run the time steps for  MS data (last 7 pond trials)
		//##################################################

		if (!runOnlyFL){

			for (int i=1; i<8; i++) {	//loop through each MS pond trial  

				for (int j=0; j<numCrabs; j++) {

					BlueCrab crab = crabs.get(i*numCrabs+j);
					if (crab.id != i) {
						System.out.println("Error -- id doesn't match with pond trial #");
						System.exit(1);
					}

					double[] tempDataMS = MSTempData.get(i-1);  //get temperature data
					int numMSTimeSteps =MSTimeSteps[i-1];  			//get the time steps for this trial

					//loop over just the days between stock and harvest of pond (i.e., age is unimportant)
					for (int k=0; k<numMSTimeSteps; k++) {

						//terminal molt for females after reach maturity
						if (crab.sex == 0 || (crab.sex == 1  && crab.CW < crab.sizeMature)) {

							//(1) accummulate degree days
							double degDays = tempDataMS[k] - 8.9;

							if (degDays > 0) 
								crab.degreeDays += degDays;

							//(2) check if molt, degree days > IP
							if (crab.degreeDays > crab.IP) {

								//Spline approach:
								meanGPM=1;
									//need to set upper bound for function, else will throw error
									if (crab.CW>maxXVal) meanGPM=spline.value(maxXVal);
									else meanGPM=spline.value(crab.CW);
								//some rounding errors in func.value; need to make sure never drops below
								if (meanGPM<lowBounds) meanGPM=lowBounds + 0.001;

								//LOGNORMAL GPM
								if (gpmErrType==0){
									logStDev=Math.sqrt(Math.log(1+(GPM_stDev*GPM_stDev)/((meanGPM-lowBounds)*(meanGPM-lowBounds))));
									logMu = Math.log(meanGPM-lowBounds) - 0.5*(logStDev*logStDev); 
									gpm = lowBounds+Math.exp(normal.nextDouble(logMu,logStDev));
								}

								//BETA GPM
								else if (gpmErrType==1) {
									//Beta distribution calculation of alpha and beta from sample mean & var
									adjMean=(meanGPM-lowBounds)/(upBounds-lowBounds);
									cv=GPM_stDev/meanGPM;	//set up the CV
									var=(cv*adjMean)*(cv*adjMean);	//calculate variance
									alpha =((1-adjMean)/var-1/adjMean)*(adjMean*adjMean);
									beta = alpha*(1/adjMean-1);								
									gpm=(upBounds-lowBounds)*betaDist.nextDouble(alpha, beta)+lowBounds;
								}
								else if (gpmErrType==2){
									gpm=normal.nextDouble(meanGPM, GPM_stDev);
								}

								crab.CW *= gpm;

								crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
										*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
										+ gamma_a*Math.pow(gamma_b,crab.CW);

								while (crab.IP > gamma_hi *  beta_a*Math.pow(beta_b,crab.CW)	) {
									crab.IP = -Math.log(uniform.nextDoubleFromTo(0, 1))
											*(beta_a*Math.pow(beta_b,crab.CW) - gamma_a*Math.pow(gamma_b,crab.CW))
											+ gamma_a*Math.pow(gamma_b,crab.CW);
								}

								crab.degreeDays = 0;  // reset to 0
							} // end check if crabs molt

						}//check if male or immature female in order to grow


						crab.ageInDays++;


					}// loop over time steps

					if (crab.sex == 0)
						mMSEst[i-1][j] = crab.CW;
					else 
						fMSEst[i-1][j-(int) (numCrabs/2.)] = crab.CW;	//females are in 2nd half of array

				}// loop over crabs
			} //loop over trials/ponds

		}//end if run only FL

		//recylce the crabs for next optimization to use for memory management
		recycleCrabs(crabs);




		//##################################################
		//calculate fitness
		//##################################################
		fitness=0;
		//compute the objective function
		double estMu=0, estStDev=0;
		double err=0;



		//########################
		//FLORIDA Males
		//########################

		//Males, sum up log likelihoods
		int lastIndex=-1;
		for (int i=0; i<mObs.length; i++) {

			int day = (int) mObs[i][0];	//get the day of observation

			int index = Arrays.binarySearch(mObsDays, day);
			if (index < 0) {
				System.out.println("Day doesn't match up, some error...");
				System.exit(0);
			}

			//only compute the mean and stdev of the predicted if it's a new day
			if (index != lastIndex){
				estMu = StatUtils.mean(mEst[index]);
				estStDev = Math.sqrt(StatUtils.variance(mEst[index]));
				lastIndex=index;
			}

			//USE full log-likelihood (normal distribution), weighted by nObservations
			//uses full LL so that the observed and predicted variances match up, thus
			//fitting the variability parameters (GPM st dev, gamma a/b) appropriately
			err+=0.5*Math.log(2.*Math.PI)+Math.log(estStDev)+0.5*((mObs[i][1]-estMu)*(mObs[i][1]-estMu))/(estStDev*estStDev);
		}
		//weight the SSE by the number of observations so equal
		//weight to males, females, and each of Darcie's pond data
		fitness+= err;



		//########################
		//Females
		//########################
		err=0;
		lastIndex=-1;

		//Females, sum up log likelihoods
		for (int i=0; i< fObs.length; i++) {
			int day = (int) fObs[i][0];	//get the day of observation
			int index = Arrays.binarySearch(fObsDays, day);
			if (index < 0) {
				System.out.println("Day doesn't match up, some error...");
				System.exit(0);
			}

			//only compute the mean and stdev of the predicted if it's a new day
			if (index != lastIndex){
				estMu = StatUtils.mean(fEst[index]);
				estStDev = Math.sqrt(StatUtils.variance(fEst[index]));
				lastIndex=index;
			}

			//USE full log-likelihood for normal distribution, weighted by nObservations
			err+=0.5*Math.log(2.*Math.PI)+Math.log(estStDev)+0.5*((fObs[i][1]-estMu)*(fObs[i][1]-estMu))/(estStDev*estStDev);
		}
		fitness+= err;



		//########################
		//Darcie's pond data
		//########################

		if (!runOnlyFL){

			// MALES
			for (int i=1; i<8; i++){ //loop over the trials (note, in the data indexing starts at trial#1)
				err=0;

				//get mean and stdev of the predicted sizes for likelihood function
				estMu = StatUtils.mean(mMSEst[i-1]);
				estStDev = Math.sqrt(StatUtils.variance(mMSEst[i-1]));

				//loop over all observations for a sex and trial #
				String key="M"+i;
				ArrayList<Double> sizes = MSObsSize.get(key);

				for (int k=0;k<sizes.size();k++) {
					//USE full log-likelihood for normal distribution, weighted by nObservations
					err+=0.5*Math.log(2.*Math.PI)+Math.log(estStDev)+0.5*((sizes.get(k)-estMu)*(sizes.get(k)-estMu))/(estStDev*estStDev);
				}
				fitness+= err;
			}


			// FEMALES
			for (int i=1; i<8; i++){ //loop over the trials (note, in the data indexing starts at trial#1)
				err=0;

				//get mean and stdev of the predicted sizes for likelihood function
				estMu = StatUtils.mean(fMSEst[i-1]);
				estStDev = Math.sqrt(StatUtils.variance(fMSEst[i-1]));

				//loop over all observations for a sex and trial #
				String key="F"+i;
				ArrayList<Double> sizes = MSObsSize.get(key);

				for (int k=0;k<sizes.size();k++) {
					//USE full log-likelihood for normal distribution, weighted by nObservations
					err+=0.5*Math.log(2.*Math.PI)+Math.log(estStDev)+0.5*((sizes.get(k)-estMu)*(sizes.get(k)-estMu))/(estStDev*estStDev);
				}
				fitness+= err;
			}
		}//end if run only FL

		if (fitness < 1){
			System.out.println("Getting negative fitness, so need a better scheme");
			System.exit(1);
		}
		fitness=1/(1+Math.log(fitness));	//want to maximize the log likelihood (how ECJ works), so back-convert from negLL
		//NOTE: can't use Koza fitness because negLL can be negative, so flip sign to maximize

		//double offset=1000000;
		//if (fitness > offset){
		//	System.out.println("Need to increase the offset for neg LL");
		//	System.exit(1);
		//}
		//fitness = Math.log(offset-fitness);	//reverse the direction and add offset so log(>0)
		//fitness = offset-fitness;	//reverse the direction and add offset so log(>0)

		((SimpleFitness)ind2.fitness).setFitness(state,
				/// ...the fitness...
				(float) fitness, 
				///... our definition of the ideal individual
				fitness == idealFitness);

		ind2.evaluated = true;
	}



	//Make these static since program will spawn multiple instances of this class for evaluate() method
	//snyrhonize to a static lock, else will use the Class as a lock and come to standstill
	protected  static ArrayList<BlueCrab> getCrabs() {
		synchronized(LOCK_1) {
			if (recycledCrabs == null || recycledCrabs.isEmpty()) return new ArrayList<BlueCrab>();
			else return recycledCrabs.remove(recycledCrabs.size()-1); //take from end so no resorting; doesn't matter with small number probably
		}
	}


	protected static void recycleCrabs(ArrayList<BlueCrab> crabs) {
		synchronized(LOCK_1) {
			if (recycledCrabs == null) recycledCrabs = new ArrayList<ArrayList<BlueCrab>>(); 
			recycledCrabs.add(crabs); //take from end so no resorting; doesn't matter with small number probably
		}
	}


}
