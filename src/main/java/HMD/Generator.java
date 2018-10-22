package HMD;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.text.DecimalFormat;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.MultimapBuilder;

public class Generator {
	public static SaturationChecker saturation;
	public static HashSet<String> uniquecheck=new HashSet<String>();
	public static List<IAtomContainer> atomextlist= new ArrayList<IAtomContainer>();
	public static List<IAtomContainer> atomsatlist= new ArrayList<IAtomContainer>();
	public static boolean verbose = false;
	static String filedir = null;
	static String molinfo= null;
	
	public static Map<String, Integer> valences; 
	static {
		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
			
		valences.put("C", 4);
		valences.put("N", 5);
		valences.put("O", 2);
		valences.put("S", 6);
		valences.put("P", 5);
		valences.put("F", 1);
		valences.put("I", 7);
		valences.put("Cl", 5);
		valences.put("Br", 5);
		valences.put("H", 1);
	}
	
	/**
	 * These are the basic sub functions used in main ones. 
	 */
	
	//This function takes a string of atom-implicit hydrogen information to build an atomcontainer
	public static IAtomContainer build(String mol) {
		IAtomContainer atomcontainer = new org.openscience.cdk.silent.AtomContainer();
		List<String> symbols = new ArrayList<String>();
        List<Integer> hydrogens = new ArrayList<Integer>();
        String[] atoms = mol.split("(?=[A-Z])");
        for (String atom : atoms) {
            String[] info = atom.split("(?=[0-9])", 2);   
            symbols.add(info[0]);
            hydrogens.add(info.length > 1 ? Integer.parseInt(info[1]):0);
        }
        
        for(int i=0;i<symbols.size();i++) {
        	atomcontainer.addAtom(new Atom(symbols.get(i)));
        	atomcontainer.getAtom(i).setImplicitHydrogenCount(hydrogens.get(i));
        }
    	return atomcontainer;
    }
	
	//Calculates the CDK canon symmetry array representing the symmetry class distribution of the molecule.
	public static long[] canonsym(IAtomContainer mol){
		int[][]	g = GraphUtil.toAdjList(mol);
	    long[] sym= Canon.symmetry(mol, g);
	    return sym;
	}
	
	//Summation of the connected bond orders.
	public static int ordsum(IAtomContainer mol, int i){
		int count=0;
		for(IBond bond: mol.getConnectedBondsList(mol.getAtom(i))){
			count=count+bond.getOrder().numeric();
	    }
		return count;
	}
	
	//Saturation checker, checking the maximum number of connected bonds of atoms.
	public static boolean satcheck(IAtomContainer mol, int i) throws CloneNotSupportedException, CDKException, IOException{
		if ((mol.getAtom(i).getImplicitHydrogenCount()+ordsum(mol,i))>= (int)valences.get(mol.getAtom(i).getSymbol())){ 
			return false;
		}else{
			return true;
		}
	}
		
	// Counting open sites of atoms.
	public static int opencounter(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		int open = valences.get(mol.getAtom(i).getSymbol()).intValue()- ordsum(mol,i) - mol.getAtom(i).getImplicitHydrogenCount(); 
		return open;
	}
	
	//It generates the InChIs of molecules.
	public static String inchigen(IAtomContainer container) throws CDKException {
		String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(container).getInchi();	
		return inchi;
	}
		
	public static final Comparator<String> ASC_ORDER = new Comparator<String>() {
	    public int compare(String e1, String e2) { 
	        return e2.compareTo(e1);
	    }
	};
	
	//The equivalent classes of molecules are ordered and enumerated in ascending order based on their open values and implicit hydrogens; as described in the paper. 
	public static ListMultimap<String,Integer> ecenumlist(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		ListMultimap<String,Integer> classes = MultimapBuilder.treeKeys(ASC_ORDER).arrayListValues().build();
		long[] sym=canonsym(acontainer);
		for(int i=0; i<acontainer.getAtomCount();i++){
			if(satcheck(acontainer, i)==true){	
				classes.put(acontainer.getAtom(i).getSymbol()+opencounter(acontainer, i)+Long.valueOf(sym[i]).intValue(), i); //The open sites and the sym values are used for labelling. Less interactions for an atom means lower sym values.
			}
		}		
		return classes;
	}
	
	// Molecule depiction generator
	public static void depict(IAtomContainer mol, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depict = new DepictionGenerator();
		depict.withCarbonSymbols().withSize(1000, 1000).withZoom(4).depict(mol).writeTo(path);
	}
	
	//clean all the IDs
	public static IAtomContainer IDclean(IAtomContainer mol){
		for(IBond bond: mol.bonds()){
			if(bond.getID()=="last"){
				bond.setID(null);
			}
			if(bond.getID()=="increased") {
				bond.setID(null);
			}
		}
	return mol;
	}
	
	//Remove the last added bond. If the order was increased, it decreases the order of the latest.
	public static void removelast(IAtomContainer mol) {
		for(IBond bnd:mol.bonds()) {
			if(bnd.getID()=="last") {
				mol.removeBond(bnd);
			}
			if(bnd.getID()=="increased") {
				BondManipulator.decreaseBondOrder(bnd);
			}
		}
	}
	
	public static List<Integer> ecindices(IAtomContainer mol) throws CloneNotSupportedException, CDKException, IOException {
		List<Integer> e =new ArrayList<Integer>();
        ListMultimap<String, Integer> ec=ecenumlist(mol);
        Object[] array=ec.keySet().toArray();
        for(int i=0;i<array.length;i++) {
        	for(int j:ec.get((String)array[i])){
        		e.add(j);
        	}
        }
        return e;
	}
	
	/**
	 * These are the main functions used for the structure generation.
	 * 
	 */
	

	/**
	 * Function is for the initialisation of the inputs and recording the duration time.
	 */
	public static void HMD(String molinfo, String filedir) throws CloneNotSupportedException, CDKException, IOException {
		long startTime = System.nanoTime(); //Recording the duration time.
		SDFWriter outFile = new SDFWriter(new FileWriter(filedir+"output.sdf"));
		List<IAtomContainer> mols= new ArrayList<IAtomContainer>();
		IAtomContainer mol=build(molinfo);
		if(verbose) {
			System.out.println("Input molecule is built and its image is stored in the given directory.");
			//depict(mol,filedir+"inputmolecule.png");
		}
        mols.add(mol);
        if(verbose) System.out.println("Start generating structures ...");
        genall(mols,ecindices(mol),outFile);
        long endTime = System.nanoTime()- startTime;
        double seconds = (double) endTime / 1000000000.0;
		DecimalFormat d = new DecimalFormat(".###");
        if(verbose) {
        	System.out.println("Number of generated structures:"+" "+uniquecheck.size());
        	System.out.println("Duration:"+" "+d.format(seconds)); //Format is second
        }
        outFile.close();	
	}
	
	/**
	 * It is the main structure Gen function saturating all the atoms of the molecule and considering all the possible extensions.
	 */
	
	public static void genall(List<IAtomContainer> mol,List<Integer> indices, SDFWriter outFile) throws CloneNotSupportedException, CDKException, IOException {
		Iterator<Integer> iterator = indices.iterator();
		while(iterator.hasNext()) {
			int index=iterator.next();
			//List<IAtomContainer> newmol = new ArrayList<IAtomContainer>(mol);
			for(IAtomContainer ml:mol) {
				atomsat(ml,index,outFile);
			}
			
			mol.addAll(atomsatlist);
			atomsatlist.clear();
			iterator.remove();
		}
	}
	public static int counts=0;
	/**
	 * This function extends the molecule until the chosen index, atom, is saturated.
	 */
	public static List<IAtomContainer> atomsat(IAtomContainer mol,int index, SDFWriter outFile) throws CloneNotSupportedException, CDKException, IOException {
		saturation = new SaturationChecker();
		List<IAtomContainer> de=atomext(mol,index);
        List<IAtomContainer> copy= new ArrayList<IAtomContainer>(de);
        atomextlist.clear();
        for(IAtomContainer ac:copy) {
        	counts++;
        	if(satcheck(ac,index)) {
        		atomsat(ac,index,outFile);
        	}else if(!satcheck(ac,index)) {
        		//saturation.isSaturated(ac) && 
        		atomsatlist.add(ac);
        		if(ConnectivityChecker.partitionIntoMolecules(ac).getAtomContainerCount() == 1 && !uniquecheck.contains(inchigen(ac))) {
        			uniquecheck.add(inchigen(ac));
        			//depict(ac,"C:\\Users\\mehme\\Desktop\\No-Backup Zone\\parallel\\"+counts+".png");
        			outFile.write(ac);
        		}
        	}
        }
        return atomsatlist;
	}
	
	/**
	 * This functions detects the target atom to add a bond between the chosen index and
	 * the others.
	 */
	public static int targetatom(ListMultimap<String, Integer> ec, String key, int index) {
		int target=0;
		List<Integer> indices=ec.get(key);
		if(indices.contains(index) && indices.size()>1) { //If size is 1 no need to consider
			if(indices.indexOf(index)!=indices.size()-1) {
				target+= indices.get(indices.indexOf(index)+1);
			}else if(indices.indexOf(index)==indices.size()-1) {
				target+= indices.get(indices.indexOf(index)-1);
			}
		}else {
			target+= indices.get(0);
		}
		return target;
	}
	
	/**
	 * The function add a bond between two atoms or increase the order of the bond.
	 */
	public static void bondadder(IAtomContainer mol, int index, int target)throws CloneNotSupportedException, CDKException, IOException {
		IBond add = mol.getBond(mol.getAtom(index), mol.getAtom(target)); 
		if(add == null){ 					
			mol.addBond(index, target, IBond.Order.SINGLE);
			mol.getBond(mol.getAtom(index), mol.getAtom(target)).setID("last");
		}
		else{
			BondManipulator.increaseBondOrder(add); 
			mol.getBond(mol.getAtom(index), mol.getAtom(target)).setID("increased"); // 
		}
	}
	
	/**
	 * This function extends the atom in the atomcontainer by adding new bond between the atom and the others.
	 */
	public static  List<IAtomContainer> atomext(IAtomContainer mol, int index) throws CloneNotSupportedException, CDKException, IOException { 	
		ListMultimap<String, Integer> ec=ecenumlist(mol);
		for(String key:ec.keySet()) {
			int target=targetatom(ec,key,index);
			if(index!=target && satcheck(mol,index) && satcheck(mol,target)){ 
				IDclean(mol); //IDs are cleaned to mark the latest bond.
				bondadder(mol,index,target);
				IAtomContainer mol2=mol.clone();
				atomextlist.add(mol2);
				removelast(mol);
			}
		}
		return atomextlist;
	}
	
	private void parseArgs(String[] args) throws ParseException
	{
		Options options = setupOptions(args);	
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			Generator.molinfo = cmd.getOptionValue("molecularinfo");
			Generator.filedir = cmd.getOptionValue("filedir");
			
			if (cmd.hasOption("verbose")) Generator.verbose = true;
		
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			HelpFormatter formatter = new HelpFormatter();
			formatter.setOptionComparator(null);
			String header = "\nGenerates structures for a given molecular information."
					+ " The input is the string of atom symbols with their number of implicit hydrogen."
					+ "For example 'C3C3C3' means three carbon atoms each of which has three implicit hydrogens."
					+ "Besides this molecular information, the directory is needed to be specified for the output"
					+ "file. \n\n";
			String footer = "\nPlease report issues at https://github.com/MehmetAzizYirik/HMD";
			formatter.printHelp( "java -jar maygen.jar", header, options, footer, true );
			throw new ParseException("Problem parsing command line");
		}
	}
	
	private Options setupOptions(String[] args)
	{
		Options options = new Options();
		Option molinfo = Option.builder("i")
			     .required(true)
			     .hasArg()
			     .longOpt("molecularinfo")
			     .desc("String of atoms with their implicit hydrogen information (required)")
			     .build();
		options.addOption(molinfo);
		Option verbose = Option.builder("v")
			     .required(false)
			     .longOpt("verbose")
			     .desc("Print messages about the duration time of the Gen")
			     .build();
		options.addOption(verbose);	
		Option filedir = Option.builder("d")
			     .required(true)
			     .hasArg()
			     .longOpt("filedir")
			     .desc("Creates and store the output sdf file in the directory (required)")
			     .build();
		options.addOption(filedir);
		return options;
	}
	
	public static void main(String[] args) throws CloneNotSupportedException, CDKException, IOException  {
		// TODO Auto-generated method stub
		Generator gen = null;
		//String[] args1= {"-i","C3C3CC2CC","-v","-d","C:\\Users\\mehme\\Desktop\\MVN\\"};
		try {
			gen = new Generator();
			gen.parseArgs(args);
			Generator.HMD(Generator.molinfo, Generator.filedir);
		} catch (Exception e) {
			// We don't do anything here. Apache CLI will print a usage text.
			if (Generator.verbose) e.getCause(); 
		}

	}
}


