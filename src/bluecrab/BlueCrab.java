package bluecrab;

import java.util.Calendar;
import java.util.TimeZone;

public class BlueCrab {

    public int sex;		    		//0 for males, 1 for females
    public int ageInDays;	    // age in days
    public double CW;		    //carapace witch
    public int shellStatus;			//1=soft, 2=peeler, 3=hard
    public int daysSinceMolt; 		//keep track of num days since last molt
    public double wt;			    //wt in grams
    public double sizeMature;	    //0=immature, 1=mature
    public double degreeDays; 	//accumulated number of degree days
    public double IP;			    //next intermolt period
    public int id;
    public int numMolts;
    public int moltsAfterMature;
    public Calendar birthday;		// birhday in Calendar representation
    public Calendar currentDate;		// current date: birthday + ageInDays
    public double popSize;				//number of crabs in this crab if it's a superindividual
    public boolean isAlive;
    public int spawnDay;
    public boolean hasSpawned;
    public double spawners;
    
    public void initialize() {
    	ageInDays = 30; 			//assume has settled out after larval period
    	degreeDays = 0;
    	numMolts=0;
    	birthday= Calendar.getInstance(TimeZone.getTimeZone("GMT"));
    	currentDate= Calendar.getInstance(TimeZone.getTimeZone("GMT"));
    	isAlive=true;
    	hasSpawned=false;
    	spawnDay=999;
    	spawners=0;
    }

    


}
