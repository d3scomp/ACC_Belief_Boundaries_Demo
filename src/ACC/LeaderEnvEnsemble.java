package ACC;

import cz.cuni.mff.d3s.deeco.annotations.KnowledgeExchange;
import cz.cuni.mff.d3s.deeco.annotations.Membership;
import cz.cuni.mff.d3s.deeco.annotations.In;
import cz.cuni.mff.d3s.deeco.annotations.Out;
import cz.cuni.mff.d3s.deeco.annotations.PeriodicScheduling;
import cz.cuni.mff.d3s.deeco.ensemble.Ensemble;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;



public class LeaderEnvEnsemble extends Ensemble {

@Membership
public static boolean membership(
		@In("coord.lGas") Double lGas,
		@In("coord.lBrake") Double lBrake,
		@In("coord.lPos") Double lPos,
		@In("coord.lSpeed") Double lSpeed,
		@In("coord.lCreationTime") Double lCreationTime,
		
		@In("member.eLGas") Double eLGas,
		@In("member.eLBrake") Double eLBrake,
		@In("member.eLPos") Double eLPos,
		@In("member.eLSpeed") Double eLSpeed,
		@In("member.eLastTime") Double eLastTime
	){
	return true;
	}

@KnowledgeExchange
@PeriodicScheduling(50)
public static void map(
		@In("coord.lGas") Double lGas,
		@In("coord.lBrake") Double lBrake,
		@Out("coord.lPos") OutWrapper<Double> lPos,
		@Out("coord.lSpeed") OutWrapper<Double> lSpeed,
		@Out("coord.lCreationTime") OutWrapper<Double> lCreationTime,
		
		@Out("member.eLGas") OutWrapper<Double> eLGas,
		@Out("member.eLBrake") OutWrapper<Double> eLBrake,
		@In("member.eLPos") Double eLPos,
		@In("member.eLSpeed") Double eLSpeed,
		@In("member.eLastTime") Double eLastTime
	) {

	eLGas.value = lGas;
	eLBrake.value = lBrake;
	lPos.value = eLPos;
	lSpeed.value = eLSpeed;
	lCreationTime.value = eLastTime;
}

}

