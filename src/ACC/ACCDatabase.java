package ACC;

import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

public class ACCDatabase {
	
	protected static final ArrayList<Double> driverSpeed = new ArrayList<Double>();
	protected static final ArrayList<Double> routeSlops = new ArrayList<Double>();
	protected static final ArrayList<Double> fTorques = new ArrayList<Double>();
	protected static final ArrayList<Double> lTorques = new ArrayList<Double>();
	protected static final ArrayList<Double> positionSeries = new ArrayList<Double>();
	protected static final ArrayList<Double> speedSeries = new ArrayList<Double>();
	protected static final double lMass = 1000;
	protected static final double fMass = 1000;
	protected static final double g = 9.80665;

	static{
		
		positionSeries.add(0.0);
		positionSeries.add(1000.0);
		positionSeries.add(2000.0);
		positionSeries.add(3000.0);
		positionSeries.add(4000.0);
		positionSeries.add(5000.0);
		positionSeries.add(6000.0);
		positionSeries.add(7000.0);
		positionSeries.add(8000.0);
		positionSeries.add(9000.0);
		positionSeries.add(10000.0);
		positionSeries.add(100000.0);

		
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(150.0);
		driverSpeed.add(170.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);
		driverSpeed.add(90.0);

		routeSlops.add(0.0);
		routeSlops.add(0.0);
		routeSlops.add(Math.PI/60);
		routeSlops.add(Math.PI/30);
		routeSlops.add(Math.PI/20);
		routeSlops.add(Math.PI/15);
		routeSlops.add(0.0);
		routeSlops.add(0.0);
		routeSlops.add(-Math.PI/18);
		routeSlops.add(-Math.PI/36);
		routeSlops.add(0.0);
		routeSlops.add(0.0);

		
		
		speedSeries.add(0.0);
		speedSeries.add(8.0);
		speedSeries.add(20.0);
		speedSeries.add(28.0);
		speedSeries.add(40.0);
		speedSeries.add(60.0);
		speedSeries.add(80.0);
		speedSeries.add(100.0);
		speedSeries.add(120.0);
		speedSeries.add(140.0);
		speedSeries.add(160.0);
		speedSeries.add(180.0);
		speedSeries.add(200.0);
		speedSeries.add(100000.0);

		
		lTorques.add(165.0);
		lTorques.add(180.0);
		lTorques.add(180.0);
		lTorques.add(170.0);
		lTorques.add(170.0);
		lTorques.add(150.0);
		lTorques.add(115.0);
		lTorques.add(97.0);
		lTorques.add(80.0);
		lTorques.add(70.0);
		lTorques.add(60.0);
		lTorques.add(50.0);
		lTorques.add(40.0);
		lTorques.add(1.0);

		fTorques.addAll(lTorques);
	}
	
	public static Double getAcceleration(Double speed, Double pos, ArrayList<Double> torques,Double gas, Double brake, Double mass){
		double FEng = gas * getValue(speedSeries, torques, speed) / 0.005;
		double FResistance = brake * 10000;
		double FEngResistance = 0.0005 * speed;
		double FHill = Math.sin(getValue(positionSeries, routeSlops, pos)) * g * mass;
		double FFinal = FEng - FResistance - FEngResistance - FHill;
		double Acceleration = FFinal / mass;
		return Acceleration;
	}
	
	
	public static Double getValue(ArrayList<Double> x, ArrayList<Double> y, Double key){
		double a[] = ArrayUtils.toPrimitive(x.toArray(new Double[x.size()]));
		double b[] = ArrayUtils.toPrimitive(y.toArray(new Double[y.size()]));
		UnivariateInterpolator interpolator = new LinearInterpolator();//Spline interpolation more accurate
		UnivariateFunction function = interpolator.interpolate(a,b);
		if(key < 0.0) return 0.0;
		double value = function.value(key);
		return value;
	}
}
