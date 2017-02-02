package japsa.bio.alignment.ppfsm;

import japsa.bio.alignment.ProbFSM.MachineState;

public class Emission {
	//This pointer is for implementing dynamic programming using a linked-list
	Emission next = null;
	//These are for backward/forward
	//public Emission fwdEmission = null;
	public Emission bwdEmission = null;


	public int gPos, mPos;
	MachineState toState;

	public double myCost;

	//String hashKey;

	public Emission(MachineState state, int mP, int gP) {
		toState = state;
		mPos = mP;
		gPos = gP;
		myCost = Double.MAX_VALUE;
	}

	static public String hashKey(String name, int mPos, int gPos) {
		return name + "_" + mPos + "_" + gPos;

	}

	public double getScore() {
		return myCost;
	}
}