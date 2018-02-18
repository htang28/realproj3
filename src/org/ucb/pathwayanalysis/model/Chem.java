package org.ucb.pathwayanalysis.model;

/**
 *
 * @author J. Christopher Anderson
 */
public class Chem {
    private final String name;
    private final String inchi;
    
    public Chem(String name, String inchi) {
        this.name = name;
        this.inchi = inchi;
    }

    public String getName() {
        return name;
    }

    public String getInchi() {
        return inchi;
    }
}
