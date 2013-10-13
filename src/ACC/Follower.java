package ACC;

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
import cz.cuni.mff.d3s.deeco.knowledge.Component;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;

public class Follower extends Component {

	public String name;
	public Double fPos = 0.0;
	public Double fSpeed = 0.0;
	public Double fGas = 0.0;
	public Double fBrake = 0.0;
	public Double fLastTime = 1.0;

	public Double fLPos = 0.0;
	public Double fLSpeed = 0.0;
	public Double fLPosMin = 0.0;
	public Double fLSpeedMin = 0.0;
	public Double fLPosMax = 0.0;
	public Double fLSpeedMax = 0.0;
	public Double fLCreationTime = 0.0;
	public Double fLTargetPos = 0.0;
	public Double fLTargetSpeed = 0.0;

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
	protected static final double THRESHOLD = 30;

	public Follower() {
		name = "F";
	}

	@Process
	@PeriodicScheduling((int) TIMEPERIOD)
	public static void computeTarget(
			@In("fPos") Double fPos,
			@In("fSpeed") Double fSpeed, 
			@In("fLPos") Double fLPos, 
			@In("fLSpeed") Double fLSpeed,
			@In("fHeadwayDistance") Double fHeadwayDistance,
			@In("fLCreationTime") Double fLCreationTime,

			@Out("fGas") OutWrapper<Double> fGas,
			@Out("fBrake") OutWrapper<Double> fBrake,

			@InOut("fLTargetPos") OutWrapper<Double> fLTargetPos,
			@InOut("fLTargetSpeed") OutWrapper<Double> fLTargetSpeed,

			@InOut("fLPosMin") OutWrapper<Double> fLPosMin,
			@InOut("fLSpeedMin") OutWrapper<Double> fLSpeedMin,
			@InOut("fLPosMax") OutWrapper<Double> fLPosMax,
			@InOut("fLSpeedMax") OutWrapper<Double> fLSpeedMax,
			@InOut("fLastTime") OutWrapper<Double> fLastTime,
			@InOut("fIntegratorError") OutWrapper<Double> fIntegratorError,
			@InOut("fErrorWindup") OutWrapper<Double> fErrorWindup
			) {
	
		computeBeliefBoundaries(fLPos, fLSpeed, fLTargetPos.value, fLCreationTime, fLPosMin, fLSpeedMin, fLPosMax, fLSpeedMax, fLastTime);
		double inaccuracy = -1;
		if (fLTargetPos.value != 0.0)
			inaccuracy = Math.max( fLPos - fLPosMin.value , fLPosMax.value - fLPos ); 

		if (inaccuracy <= THRESHOLD) {
			fLTargetPos.value = fLPos;
			fLTargetSpeed.value = fLSpeed;
			computeTargetByCACC(fPos, fSpeed, fLTargetPos.value, fLTargetSpeed.value, fGas, fBrake, fIntegratorError, fErrorWindup);
		} else {
			if ((fLPos - fPos) <= fHeadwayDistance) {
				fLTargetPos.value = fLPos;
				fLTargetSpeed.value = fLSpeed;
			} else {
				System.out.println("ACC _____ no leader.");
				fLTargetPos.value = fPos + fHeadwayDistance;
				fLTargetSpeed.value = DESIRED_SPEED;
			}
			computeTargetByACC(fLPos, fPos, fLTargetPos.value, fGas, fBrake);
		}
	}
	

	
	private static void computeBeliefBoundaries(
			@In("fLPos") Double fLPos,
			@In("fLSpeed") Double fLSpeed,
			@In("fLTargetPos") Double fLTargetPos,
			@In("fLCreationTime") Double fLCreationTime,

			@InOut("fLPosMin") OutWrapper<Double> fLPosMin,
			@InOut("fLSpeedMin") OutWrapper<Double> fLSpeedMin,
			@InOut("fLPosMax") OutWrapper<Double> fLPosMax,
			@InOut("fLSpeedMax") OutWrapper<Double> fLSpeedMax,
			@InOut("fLastTime") OutWrapper<Double> fLastTime) {

		double currentTime = System.nanoTime() / SEC_NANOSEC_FACTOR;
		double[] minBoundaries = new double[1];
		double[] maxBoundaries = new double[1];
		double startTime = 0.0;

		if(fLTargetPos != 0.0 ) {

			if (fLCreationTime <= fLastTime.value) {
				startTime = fLastTime.value;
			} else {
				startTime = fLCreationTime;
				fLPosMin.value = fLPos;
				fLPosMax.value = fLPos;
				fLSpeedMin.value = fLSpeed;
				fLSpeedMax.value = fLSpeed;
			}

			// ---------------------- knowledge evaluation --------------------------------

			double accMin = ACCDatabase.getAcceleration(fLSpeedMin.value,
					fLPosMin.value, ACCDatabase.lTorques, 0.0, 1.0,
					ACCDatabase.lMass);
			double accMax = ACCDatabase.getAcceleration(fLSpeedMax.value,
					fLPosMax.value, ACCDatabase.lTorques, 1.0, 0.0,
					ACCDatabase.lMass);

			FirstOrderIntegrator integrator = new MidpointIntegrator(1);
			integrator.setMaxEvaluations((int) TIMEPERIOD);
			FirstOrderDifferentialEquations f = new Derivation(); 
			// ------------- min ----------------------

			minBoundaries[0] = accMin;
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLSpeedMin.value += minBoundaries[0];
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLPosMin.value += minBoundaries[0];
			// ------------- max ----------------------

			maxBoundaries[0] = accMax;
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLSpeedMax.value += maxBoundaries[0];
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLPosMax.value += maxBoundaries[0];

			System.out.println("//... pos: min " + fLPosMin.value + " ... max "	+ fLPosMax.value + " time :" + currentTime);
			System.out.println("//... speed: min " + fLSpeedMin.value + " ... max " + fLSpeedMax.value);
		}
		fLastTime.value = currentTime;
	}


	private static void computeTargetByCACC(
		@In("fPos") Double fPos,
		@In("fSpeed") Double fSpeed, 
		@In("fLTargetPos") Double fLTargetPos,
		@In("fLTargetSpeed") Double fLTargetSpeed,

		@Out("fGas") OutWrapper<Double> fGas,
		@Out("fBrake") OutWrapper<Double> fBrake,

		@InOut("fIntegratorError") OutWrapper<Double> fIntegratorError,
		@InOut("fErrorWindup") OutWrapper<Double> fErrorWindup) {

		System.out.println("CACC ____ takes the pos and the speed from wirless connection.");
		if (fLTargetPos == 0.0) {
			
			fGas.value = 0.0;
			fBrake.value = 0.0;
			
		} else {
			
			double timePeriodInSeconds = TIMEPERIOD / SEC_MILISEC_FACTOR;
			double distanceError = -DESIRED_DISTANCE + fLTargetPos - fPos;
			double pidDistance = KP_D * distanceError;
			double error = pidDistance + fLTargetSpeed - fSpeed;
			fIntegratorError.value += (KI_S * error + KT_S * fErrorWindup.value) * timePeriodInSeconds;
			double pidSpeed = KP_S * error + fIntegratorError.value;
			fErrorWindup.value = saturate(pidSpeed) - pidSpeed;
	
			if (pidSpeed >= 0) {
				fGas.value = pidSpeed;
				fBrake.value = 0.0;
			} else {
				fGas.value = 0.0;
				fBrake.value = -pidSpeed;
			}
		}
	}
	

	private static void computeTargetByACC(		
			@In("fLPos") Double fLPos,
			@In("fPos") Double fPos,
			@In("fLTargetPos") Double fLTargetPos,

			@Out("fGas") OutWrapper<Double> fGas,
			@Out("fBrake") OutWrapper<Double> fBrake
			) {

			System.out.println("ACC _____ takes the pos and the speed from the headway sensors.");
			if (fLTargetPos == 0.0) {
				
				fGas.value = 0.0;
				fBrake.value = 0.0;
				
			} else {
				if ((fLPos - fPos) >= DESIRED_DISTANCE) {
					fGas.value = 1.0;
					fBrake.value = 0.0;
				} else {
					fGas.value = 0.0;
					fBrake.value = 1.0;
				}
			}
	}

	
	private static double saturate(double val) {
		if (val > 1)
			val = 1;
		else if (val < -1)
			val = -1;
		return val;
	}

	
	private static class Derivation implements FirstOrderDifferentialEquations {

		@Override
		public int getDimension() {
			return 1;
		}

		@Override
		public void computeDerivatives(double t, double[] y, double[] yDot)
				throws MaxCountExceededException, DimensionMismatchException {
			int params = 1;
			int order = 1;
			DerivativeStructure x = new DerivativeStructure(params, order, 0,
					y[0]);
			DerivativeStructure f = x.divide(t);
			yDot[0] = f.getValue();
		}
	}
}