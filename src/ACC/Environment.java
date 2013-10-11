package ACC;


import cz.cuni.mff.d3s.deeco.annotations.In;
import cz.cuni.mff.d3s.deeco.annotations.Out;
import cz.cuni.mff.d3s.deeco.annotations.InOut;
import cz.cuni.mff.d3s.deeco.annotations.PeriodicScheduling;
import cz.cuni.mff.d3s.deeco.annotations.Process;
import cz.cuni.mff.d3s.deeco.knowledge.Component;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;


public class Environment extends Component {

	public String eName = "E";
	public Double eFGas = 0.0;    
	public Double eFBrake = 0.0;  
	public Double eLGas = 0.0;      
	public Double eLBrake = 0.0;    
	
	public Double eFPos = 0.0;
	public Double eFSpeed = 0.0; 
	public Double eLPos = 60.0;
	public Double eLSpeed = 0.0;
	public Double eLastTime = 0.0;

	
	protected static final double TIMEPERIOD = 100;
	protected static final double SEC_NANOSECOND_FACTOR = 1000000000;
	protected static final double SEC_MILISEC_FACTOR = 1000;
	
	
	
	public Environment() {
	}	
	
	@Process 
	@PeriodicScheduling((int) TIMEPERIOD)
	public static void environmentResponse(
			@In("eLGas") Double eLGas,
			@In("eLBrake") Double eLBrake,
			@In("eFGas") Double eFGas,
			@In("eFBrake") Double eFBrake,
			
			@InOut("eFPos") OutWrapper<Double> eFPos,			 
			@InOut("eFSpeed") OutWrapper<Double> eFSpeed,
			@InOut("eLPos") OutWrapper<Double> eLPos,			 
			@InOut("eLSpeed") OutWrapper<Double> eLSpeed,
			
			@Out("eLastTime") OutWrapper<Double> eLastTime
			){
	
		double currentTime = System.nanoTime()/SEC_NANOSECOND_FACTOR;
		double timePeriodInSeconds = TIMEPERIOD/SEC_MILISEC_FACTOR;
		
		// ----------------------- leader ----------------------------------------------------------------------
		double lAcceleration = ACCDatabase.getAcceleration(eLSpeed.value, eLPos.value, ACCDatabase.lTorques, eLGas, eLBrake, ACCDatabase.lMass);
		eLSpeed.value += lAcceleration * timePeriodInSeconds;
		eLPos.value += eLSpeed.value * timePeriodInSeconds;
 		//------------------------ follower ---------------------------------------------------------------------
		double fAcceleration = ACCDatabase.getAcceleration(eFSpeed.value, eFPos.value, ACCDatabase.fTorques, eFGas, eFBrake,ACCDatabase.fMass);
		eFSpeed.value += fAcceleration * timePeriodInSeconds; 
		eFPos.value += eFSpeed.value * timePeriodInSeconds;
 		//--------------------------------------------------------------------------------------------------------
		
		eLastTime.value = currentTime;
		System.out.println("Speed leader : "+eLSpeed.value+", pos : "+eLPos.value+"... time :"+currentTime);
		System.out.println("Speed follower : "+eFSpeed.value+", pos : "+eFPos.value+"... time :"+currentTime);
		System.out.println("................... distance : "+(eLPos.value - eFPos.value));
	}
}