package ACC;


import cz.cuni.mff.d3s.deeco.annotations.KnowledgeExchange;
import cz.cuni.mff.d3s.deeco.annotations.Membership;
import cz.cuni.mff.d3s.deeco.annotations.In;
import cz.cuni.mff.d3s.deeco.annotations.Out;
import cz.cuni.mff.d3s.deeco.annotations.PeriodicScheduling;
import cz.cuni.mff.d3s.deeco.ensemble.Ensemble;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;



public class FollowerEnvEnsemble extends Ensemble {

	@Membership
	public static boolean membership(
			@In("coord.fPos") Double fPos,
			@In("coord.fSpeed") Double fSpeed,
			@In("coord.fGas") Double fGas,
			@In("coord.fBrake") Double fBrake,
		
			@In("member.eFGas") Double eFGas,
			@In("member.eFBrake") Double eFBrake,
			@In("member.eFPos") Double eFPos,
			@In("member.eFSpeed") Double eFSpeed
		){
			return true;
	}
	
	@KnowledgeExchange
	@PeriodicScheduling(50)
	public static void map(
			@Out("coord.fPos")  OutWrapper<Double> fPos,
			@Out("coord.fSpeed")  OutWrapper<Double> fSpeed,
			@In("coord.fGas") Double fGas,
			@In("coord.fBrake") Double fBrake,
		
			@Out("member.eFGas")  OutWrapper<Double> eFGas,
			@Out("member.eFBrake")  OutWrapper<Double> eFBrake,
			@In("member.eFPos") Double eFPos,
			@In("member.eFSpeed") Double eFSpeed
	) {
	
		eFGas.value = fGas;
		eFBrake.value = fBrake;
		fPos.value = eFPos;
		fSpeed.value = eFSpeed;
	}
}