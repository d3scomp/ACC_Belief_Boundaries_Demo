package ACC;


import java.util.ArrayList;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.MidpointIntegrator;


import cz.cuni.mff.d3s.deeco.annotations.In;
import cz.cuni.mff.d3s.deeco.annotations.InOut;
import cz.cuni.mff.d3s.deeco.annotations.Out;
import cz.cuni.mff.d3s.deeco.annotations.PeriodicScheduling;
import cz.cuni.mff.d3s.deeco.annotations.Process;
import cz.cuni.mff.d3s.deeco.annotations.TriggerOnChange;
import cz.cuni.mff.d3s.deeco.knowledge.Component;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;

public class Follower extends Component {

	public String name;
	public Double fPos = 0.0;      
	public Double fSpeed = 0.0;    
	public Double fGas = 0.0 ;
	public Double fBrake = 0.0;
	public Double fLastTime = 0.0;
	
	public Double fLPos = 0.0;       
	public Double fLSpeed = 0.0;  
	public Double fLPosMin = 0.0;       
	public Double fLSpeedMin = 0.0;     
	public Double fLPosMax = 0.0;       
	public Double fLSpeedMax = 0.0;     
	public Double fLCreationTime = 0.0;	  
	public Double fLTargetPos = 0.0;       
	public Double fLTargetSpeed = 0.0;  

	public Double fInaccuracy = -1.0;
	public Double fHeadwayDistance = 100.0;
	public Double fIntegratorError = 0.0;
	public Double fErrorWindup = 0.0;
	
	
	protected static final double KP_D = 0.193;
	protected static final double KP_S = 0.12631;
	protected static final double KI_S = 0.001;
	protected static final double KT_S = 0.01;
	protected static final double SEC_NANOSEC_FACTOR = 1000000000;
	protected static final double TIMEPERIOD = 100;
	protected static final double SEC_MILISEC_FACTOR = 1000;
	protected static final double DESIRED_DISTANCE = 50;
	protected static final double DESIRED_SPEED = 90;
	protected static final double THRESHOLD = 15;

    
	
	public Follower() {
		name = "F";
	}
	
	
	@Process
	@PeriodicScheduling((int) TIMEPERIOD)
	public static void speedControl(
				@In("fPos") Double fPos,
				@In("fSpeed") Double fSpeed,
				@In("fLTargetPos") Double fLTargetPos,
				@In("fLTargetSpeed") Double fLTargetSpeed,
				@In("fInaccuracy") Double fInaccuracy,
				
				@Out("fGas") OutWrapper<Double> fGas,
				@Out("fBrake") OutWrapper<Double> fBrake,
				
				@InOut("fIntegratorError") OutWrapper<Double> fIntegratorError,
				@InOut("fErrorWindup") OutWrapper<Double> fErrorWindup
			) {
	
			double timePeriodInSeconds = TIMEPERIOD/SEC_MILISEC_FACTOR;
			double distanceError = - DESIRED_DISTANCE + fLTargetPos - fPos;
			double pidDistance = KP_D * distanceError;
			double error = pidDistance + fLTargetSpeed - fSpeed;
			fIntegratorError.value += (KI_S * error + KT_S * fErrorWindup.value) * timePeriodInSeconds;
			double pidSpeed = KP_S * error + fIntegratorError.value;
			fErrorWindup.value = saturate(pidSpeed) - pidSpeed;

			if( fInaccuracy == -1.0){
				fGas.value = 0.0;
				fBrake.value = 0.0;
			}else{
				if(pidSpeed >= 0){
				fGas.value = pidSpeed;
				fBrake.value = 0.0;
				}else{
					fGas.value = 0.0;
					fBrake.value = -pidSpeed;
				}
			}
	}
	
	
	
	@Process
	@PeriodicScheduling((int) TIMEPERIOD)
	public static void computeTarget(
			@In("fPos") Double fPos,
			@In("fLPos") Double fLPos,
			@In("fLSpeed") Double fLSpeed,
			@In("fLPosMin") Double fLPosMin,
			@In("fLSpeedMin") Double fLSpeedMin,
			@In("fInaccuracy") Double fInaccuracy, // Trigger
			@In("fHeadwayDistance") Double fHeadwayDistance,

			@InOut("fLTargetPos") OutWrapper<Double> fLTargetPos, // InOut : to not have null if we did not enter if condition
			@InOut("fLTargetSpeed") OutWrapper<Double> fLTargetSpeed
            ){
		if ( fInaccuracy <= THRESHOLD){
			computeTargetByCACC();
			fLTargetPos.value = fLPos;
			fLTargetSpeed.value =  fLSpeed;
		} else {
			if( (fLPos - fPos) <= fHeadwayDistance ){
				computeTargetByACC();
				fLTargetPos.value = fLPos;
				fLTargetSpeed.value =  fLSpeed;
			}else{
				fLTargetPos.value = fPos + fHeadwayDistance;
				fLTargetSpeed.value =  DESIRED_SPEED;
				System.out.println("ACC   _____   no leader.");
			}
		}
 	}
	
	
	
	private static void computeTargetByCACC(){
		System.out.println("CACC ____ takes the pos and the speed from wirless connection.");
 	}
	
	
	private static void computeTargetByACC(){
		System.out.println("ACC   _____  takes the pos and the speed from the headway sensors.");
	}
	
	
	@Process
	@PeriodicScheduling((int) TIMEPERIOD)
	public static void computeBeliefBoundaries(
			@In("fLPos") Double fLPos,
			@In("fLSpeed") Double fLSpeed,
			@In("fLTargetPos")  Double fLTargetPos,		      
			@In("fLCreationTime") Double fLCreationTime,
			
			@InOut("fInaccuracy") OutWrapper<Double> fInaccuracy,

			@InOut("fLPosMin") OutWrapper<Double> fLPosMin,
			@InOut("fLSpeedMin") OutWrapper<Double> fLSpeedMin,
			@InOut("fLPosMax") OutWrapper<Double> fLPosMax,
			@InOut("fLSpeedMax") OutWrapper<Double> fLSpeedMax,
			@InOut("fLastTime") OutWrapper<Double> fLastTime
			){
		
		double currentTime = System.nanoTime()/SEC_NANOSEC_FACTOR;
		double[] minBoundaries = new double[1]; 
		double[] maxBoundaries = new double[1];
		double startTime = 0.0;
		
		if( fLCreationTime <= fLastTime.value ){
			startTime = fLastTime.value;
		}else{
			startTime = fLCreationTime;
			fLPosMin.value = fLPos;
			fLPosMax.value = fLPos;
			fLSpeedMin.value = fLSpeed;
			fLSpeedMax.value = fLSpeed;
		}
		
		//------------------------------------------------ knowledge evaluation ------------------------------------------
		try{
			
			double accMin = ACCDatabase.getAcceleration(fLSpeedMin.value, fLPosMin.value, ACCDatabase.lTorques, 0.0, 1.0, ACCDatabase.lMass);
			double accMax = ACCDatabase.getAcceleration(fLSpeedMax.value, fLPosMax.value, ACCDatabase.lTorques, 1.0, 0.0, ACCDatabase.lMass);
			
			FirstOrderIntegrator integrator = new MidpointIntegrator(1);
			integrator.setMaxEvaluations((int) TIMEPERIOD);
			FirstOrderDifferentialEquations f = new Derivation(); // why I should put if F^min and F^max
			//------------- min ----------------------

			minBoundaries[0] = accMin;
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLSpeedMin.value += minBoundaries[0];
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLPosMin.value += minBoundaries[0];
			//------------- max ----------------------
			
			maxBoundaries[0] = accMax;
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLSpeedMax.value += maxBoundaries[0];
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLPosMax.value += maxBoundaries[0];

			
			System.out.println("//... pos: min "+fLPosMin.value+" ... max "+fLPosMax.value+"      time :"+currentTime);
			System.out.println("//... speed: min "+fLSpeedMin.value+" ... max "+fLSpeedMax.value);
		
		} catch ( Exception e ){
			System.err.println("error : "+e.getMessage()); //the error at the first of the execution because of the fLastTime is zero => the integration range is so big   
		}
		
		if( fLTargetPos == 0.0 )
			fInaccuracy.value = -1.0;
		else
			fInaccuracy.value = fLPos - fLPosMin.value; // do I put the inaccuracy for both min/max, fInaccuracy = Math.max( fLPos - fLPosMin.value , fLPosMax.value - fLPos ); or only what I care about which is the min boundary?

		fLastTime.value = currentTime;
	}

	
	private static double saturate(double val) {
		if(val > 1) val = 1;
		else if(val < -1) val = -1;
		return val;
	}
	
	
	private static class Derivation implements FirstOrderDifferentialEquations{

		@Override
		public int getDimension() {
			// TODO Auto-generated method stub
			return 1;
		}

		@Override
		public void computeDerivatives(double t, double[] y, double[] yDot)
				throws MaxCountExceededException, DimensionMismatchException {
			// TODO Auto-generated method stub
			int params = 1;
			int order = 1;
			DerivativeStructure x = new DerivativeStructure(params, order, 0, y[0]);
			DerivativeStructure f = x.divide(t);
			yDot[0] = f.getValue();
		}
	}
}
