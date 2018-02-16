/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.andersonlab.terminalolefins.stoichiometry.model;

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
