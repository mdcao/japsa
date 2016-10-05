package japsa.bio.hts.scaffold;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.bio.hts.scaffold.Node;
import japsa.bio.hts.scaffold.Vertex;

public 	class Node{
	Vertex v;
	boolean dir;
	Node(Vertex v, boolean dir){
		this.v=v;
		this.dir=dir;
	}
	public Vertex getVertex(){
		return v;
	}
	public void setVertex(Vertex v){
		this.v = v;
	}
	public boolean getDirection(){
		return dir;
	}
	public void setDirection(boolean dir){
		this.dir=dir;
	}
	public Node getRC(){
		return new Node(v,!dir);
	}
	public Sequence getSeq(){
		return dir?v.getSequence():Alphabet.DNA.complement(v.getSequence());
	}
	public String toString(){
		return v.getLabel()+ (dir?"+":"-");
	}
}
