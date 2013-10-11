package ACC;

import cz.cuni.mff.d3s.deeco.annotations.KnowledgeExchange;
import cz.cuni.mff.d3s.deeco.annotations.Membership;
import cz.cuni.mff.d3s.deeco.annotations.In;
import cz.cuni.mff.d3s.deeco.annotations.Out;
import cz.cuni.mff.d3s.deeco.annotations.PeriodicScheduling;
import cz.cuni.mff.d3s.deeco.ensemble.Ensemble;
import cz.cuni.mff.d3s.deeco.knowledge.OutWrapper;



public class FollowerLeaderEnsemble extends Ensemble {

	@Membership
	public static boolean membership(
			@In("coord.fLPos") Double fLPos,
			@In("coord.fLSpeed") Double fLSpeed,
			@In("coord.fLCreationTime") Double fLCreationTime,
			@In("coord.fHeadwayDistance") Double fHeadwayDistance,
	
			@In("member.lPos") Double lPos,
			@In("member.lSpeed") Double lSpeed,
			@In("member.lCreationTime") Double lCreationTime
		){
		if( (fLPos - lPos) <= 2*fHeadwayDistance )
			return true;
		return false;
	}
	
	@KnowledgeExchange
	@PeriodicScheduling(50)
	public static void map(
			@Out("coord.fLPos") OutWrapper<Double> fLPos,
			@Out("coord.fLSpeed") OutWrapper<Double> fLSpeed,
			@Out("coord.fLCreationTime") OutWrapper<Double> fLCreationTime,
	
			@In("member.lPos") Double lPos,
			@In("member.lSpeed") Double lSpeed,
			@In("member.lCreationTime") Double lCreationTime
	) {
	
	fLPos.value = lPos;
	fLSpeed.value = lSpeed;
	fLCreationTime.value = lCreationTime;
		
	}

}


