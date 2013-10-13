package ACC;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.MidpointIntegrator;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

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

	public Double fLPos = 0.0;
	public Double fLSpeed = 0.0;
	public Double fLCreationTime = 0.0;
	public Double fLTargetPos = 0.0;
	public Double fLTargetSpeed = 0.0;
	public Double fHeadwayDistance = 100.0;

	protected static double fLastTime = 0.0;
	protected static double fLPosMin = 0.0;
	protected static double fLSpeedMin = 0.0;
	protected static double fLPosMax = 0.0;
	protected static double fLSpeedMax = 0.0;
	protected static double fIntegratorError = 0.0;
	protected static double fErrorWindup = 0.0;

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
			@InOut("fLTargetSpeed") OutWrapper<Double> fLTargetSpeed
			) {
	
		computeBeliefBoundaries(fLPos, fLSpeed, fLTargetPos.value, fLCreationTime );
		double inaccuracy = -1;
		if (fLTargetPos.value != 0.0)
			inaccuracy = Math.max( fLPos - fLPosMin, fLPosMax - fLPos ); 

		if (inaccuracy <= THRESHOLD) {
			computeTargetByCACC();
			fLTargetPos.value = fLPos;
			fLTargetSpeed.value = fLSpeed;
		} else {
			if ((fLPos - fPos) <= fHeadwayDistance) {
				computeTargetByACC();
				fLTargetPos.value = fLPos;
				fLTargetSpeed.value = fLSpeed;
			} else {
				fLTargetPos.value = fPos + fHeadwayDistance;
				fLTargetSpeed.value = DESIRED_SPEED;
				System.out.println("ACC _____ no leader.");
			}
		}
		ArrayList<Double> ar=speedControl(fPos, fSpeed, fLTargetPos.value, fLTargetSpeed.value);
		fGas.value = ar.get(0);
		fBrake.value = ar.get(1);
	}
	

	
	private static ArrayList<Double> speedControl( Double fPos, Double fSpeed, Double fLTargetPos, Double fLTargetSpeed ) {

		ArrayList<Double> result = new ArrayList<Double>();
		if (fLTargetPos == 0.0) {
			result.add(0.0);
			result.add(0.0);
		} else {
			double timePeriodInSeconds = TIMEPERIOD / SEC_MILISEC_FACTOR;
			double distanceError = -DESIRED_DISTANCE + fLTargetPos - fPos;
			double pidDistance = KP_D * distanceError;
			double error = pidDistance + fLTargetSpeed - fSpeed;
			fIntegratorError += (KI_S * error + KT_S * fErrorWindup)
					* timePeriodInSeconds;
			double pidSpeed = KP_S * error + fIntegratorError;
			fErrorWindup = saturate(pidSpeed) - pidSpeed;

			if (pidSpeed >= 0) {
				result.add(pidSpeed);
				result.add(0.0);
			} else {
				result.add(0.0);
				result.add(-pidSpeed);
			}
		}
		
		return result;
	}


	private static void computeBeliefBoundaries( Double fLPos, Double fLSpeed, Double fLTargetPos, Double fLCreationTime ) {

		double currentTime = System.nanoTime() / SEC_NANOSEC_FACTOR;
		double[] minBoundaries = new double[1];
		double[] maxBoundaries = new double[1];
		double startTime = 0.0;

		if(fLTargetPos != 0.0 ) {

			if (fLCreationTime <= fLastTime) {
				startTime = fLastTime;
			} else {
				startTime = fLCreationTime;
				fLPosMin = fLPos;
				fLPosMax = fLPos;
				fLSpeedMin = fLSpeed;
				fLSpeedMax = fLSpeed;
			}

			// ---------------------- knowledge evaluation --------------------------------

			double accMin = ACCDatabase.getAcceleration(fLSpeedMin,
					fLPosMin, ACCDatabase.lTorques, 0.0, 1.0,
					ACCDatabase.lMass);
			double accMax = ACCDatabase.getAcceleration(fLSpeedMax,
					fLPosMax, ACCDatabase.lTorques, 1.0, 0.0,
					ACCDatabase.lMass);

			FirstOrderIntegrator integrator = new MidpointIntegrator(1);
			integrator.setMaxEvaluations((int) TIMEPERIOD);
			FirstOrderDifferentialEquations f = new Derivation(); 
			// ------------- min ----------------------

			minBoundaries[0] = accMin;
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLSpeedMin += minBoundaries[0];
			integrator.integrate(f, startTime, minBoundaries, currentTime, minBoundaries);
			fLPosMin += minBoundaries[0];
			// ------------- max ----------------------

			maxBoundaries[0] = accMax;
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLSpeedMax += maxBoundaries[0];
			integrator.integrate(f, startTime, maxBoundaries, currentTime, maxBoundaries);
			fLPosMax += maxBoundaries[0];

			System.out.println("//... pos: min " + fLPosMin + " ... max "	+ fLPosMax + " time :" + currentTime);
			System.out.println("//... speed: min " + fLSpeedMin + " ... max " + fLSpeedMax);
		}
		fLastTime = currentTime;
	}


	private static void computeTargetByCACC() {
		System.out.println("CACC ____ takes the pos and the speed from wirless connection.");
	}
	

	private static void computeTargetByACC() {
		System.out.println("ACC _____ takes the pos and the speed from the headway sensors.");
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