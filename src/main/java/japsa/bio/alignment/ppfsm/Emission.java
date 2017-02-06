package japsa.bio.alignment.ppfsm;

import japsa.bio.alignment.ppfsm.state.MachineState;

public class Emission {
	//This pointer is for implementing dynamic programming using a linked-list
	Emission next = null;
	
	//These are for backward/forward
	public Emission fwdEmission = null;
	public Emission bwdEmission = null;

	public int gPos, iteration;
	MachineState toState;

	public double myCost;

	//String hashKey;

	public Emission(MachineState state, int gP, int iter) {
		toState = state;
		gPos = gP;
		iteration = iter;		
		myCost = Double.MAX_VALUE;
	}

	static public String hashKey(MachineState state, int gPos, int iter) {
		return state.getName() + "_" + gPos + "_" + iter;

	}

	public double getScore() {
		return myCost;
	}
}