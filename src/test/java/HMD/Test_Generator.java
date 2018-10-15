package HMD;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertTrue;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import static org.junit.Assert.assertEquals;

import java.io.FileWriter;
import java.io.IOException;

/**
 * @cdk.module test-standard
 */

public class Test_Generator {
	private IAtomContainer atomContainer;
	public  String Desktop=(System.getProperty("user.home") + "\\Desktop").replace("\\", "/");
	@Before
	/**
	 * This function is for clearing the lists if there are something stored
	 * from the tests of the other functions.
	 */
	public void cleanlists() {
		Generator.uniquecheck.clear();
		Generator.atomextlist.clear();
		Generator.atomsatlist.clear();
	}
	
	@Before
	
	/**
	 * This atomcontainer is built as the example input for all the test cases.
	 * @throws Exception
	 */
	
    public void setUp() throws Exception {
		IAtomContainer ac = new org.openscience.cdk.silent.AtomContainer();

		ac.addAtom(new Atom("C"));
		ac.addAtom(new Atom("C"));
		ac.addAtom(new Atom("C"));
		ac.addAtom(new Atom("C"));
		ac.addAtom(new Atom("C"));
		ac.addAtom(new Atom("C"));

		ac.getAtom(0).setImplicitHydrogenCount(3);
		ac.getAtom(1).setImplicitHydrogenCount(3);
		ac.getAtom(2).setImplicitHydrogenCount(2);
		ac.getAtom(3).setImplicitHydrogenCount(2);
		ac.getAtom(4).setImplicitHydrogenCount(1);
		ac.getAtom(5).setImplicitHydrogenCount(1);
		
		this.atomContainer=ac;
    }
	
	@Test
	/**
	 * This function simply takes a string including the atom-implicit hydrogen information.
	 * such as "C1C1C1". Thus the atom container should include three carbons. Each carbon atom
	 * is attached with a implicit hydrogen.
	 */
	public void test_build() throws CloneNotSupportedException, CDKException, IOException {
		String molecule="C3C3C2C2C1C1";
		IAtomContainer acon=Generator.build(molecule);
		//The number of atoms: 6
 		assertEquals(6,acon.getAtomCount());
		//The number of implicit hydrogens: 12
 		int count=0;
 		for(int i=0;i<acon.getAtomCount();i++) {
 			count+=acon.getAtom(i).getImplicitHydrogenCount();
 		}
 		assertEquals(12,count);
	}
	
	@Test
	/**
	 * Function calculates equivalence classes using canonsym, then classifies the atoms.
	 * The equivalence classes are enumerated based on implicit hydrogens and symmetry values
	 * of the atoms. 
	 * 
	 */
	
    public void test_ecenumlist() throws CloneNotSupportedException, CDKException, IOException {	
		ListMultimap<String, Integer> map= ArrayListMultimap.create();
		map.put("C31", 4);
		map.put("C31", 5);
		map.put("C23", 2);
		map.put("C23", 3);
		map.put("C15", 0);
		map.put("C15", 1);
		
		//System.out.println(Functions.ecenum(atomContainer));
        assertEquals(map,Generator.ecenumlist(atomContainer));
    }
	
	@Test
	/**
	 * The canon values are calculated by using CDK Canon class. That calculates the symmetry classes.
	 */
	
	public void test_canonsym() {
		long[] sym= {5,5,3,3,1,1};
		assertEquals(Arrays.toString(sym), Arrays.toString(Generator.canonsym(atomContainer)));
	}
	
	@Test
	/**
	 * This function is for checking the saturation of an atom in an atomcontainer.
	 * If the atom is not saturated, the function returns boolean value.
	 */
	public void test_satcheck() throws CloneNotSupportedException, CDKException, IOException {
		assertTrue(Generator.satcheck(atomContainer, 0)==true);
	}
	
	@Test
	/**
	 * The function calculates the open sites of an atom. Here, atom at the 1th index has 1 open site.
	 */
	public void test_opencounter() throws CloneNotSupportedException, CDKException, IOException {
		//System.out.println(Functions.opencounter(atomContainer, 1));
		/**
		 * The atom with index 1 has only one open site remained.
		 */
		assertEquals(1,Generator.opencounter(atomContainer, 1));
	}
	
	
	
	/**
	 * The function is simply for depiction of molecules.
	 *
	public void test_depict() throws CloneNotSupportedException, CDKException, IOException {
		Generator.depict(atomContainer,Desktop+"\\example.png"); 
	}**/
	
	@Test
	/**
	 * The function is for the bond extension. For the input atom container and the index,
	 * the atom container is extended from the given index by considering all the possible 
	 * target atoms from different equivalence classes. 
	 **/
	public void test_atomext() throws CloneNotSupportedException, CDKException, IOException {
		cleanlists();
		IAtomContainer mol=Generator.build("C3C3C2C2C1C1CCC");
		//By calling the atomext for the index 0, 4 structures are expected to be generated.
		assertEquals(4,Generator.atomext(mol,0).size());
	}
	
	@Test
	/**
	 * Atomsat function reimplements atomext function until the atom of the chosen index is saturated.
	 * 
	 */
	public void test_atomsat() throws CloneNotSupportedException, CDKException, IOException {
		cleanlists();
		IAtomContainer mol=Generator.build("C2C2C1C1C1C1");
		SDFWriter outFile = new SDFWriter(new FileWriter(Desktop+"\\testoutput.sdf"));
		//By calling the atomsat for the index 0, 5 structures are expected to be generated.
		assertEquals(5,Generator.atomsat(mol,0,outFile).size());
	}
	
	@Test
	/**
	 * The genall function generates all the possible structures by saturating each of the indices
	 * grouped in equivalence classes.
	 */
	public void test_genall()  throws CloneNotSupportedException, CDKException, IOException{
		cleanlists();
		SDFWriter outFile = new SDFWriter(new FileWriter(Desktop+"\\testoutput.sdf"));
		List<IAtomContainer> mols= new ArrayList<IAtomContainer>();
		IAtomContainer mol=Generator.build("C3C3C2C2C1C1");
		mols.add(mol);
		Generator.genall(mols,Generator.ecindices(mol),outFile);
		/**By calling the genall function, the number of possible extensions is 1032.
		 * Among these extensions, the saturated ones are writen in to the output file.
		 */
		assertEquals(1032,mols.size());
	}
	
}

