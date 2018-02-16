/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.andersonlab.terminalolefins.stoichiometry.model;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author J. Christopher Anderson
 */
public class Pathway {
    private final String name;
    private final  Set<Chem> inputs;
    private final  Set<Chem> outputs;
    private final  List<Rxn> reactions;
    private final  Map<String, Chem> chemicals;

    public Pathway(String name, Set<Chem> inputs, Set<Chem> outputs, List<Rxn> reactions, Map<String, Chem> chemicals) {
        this.name = name;
        this.inputs = inputs;
        this.outputs = outputs;
        this.reactions = reactions;
        this.chemicals = chemicals;
    }

    public String getName() {
        return name;
    }

    public Set<Chem> getInputs() {
        return inputs;
    }

    public Set<Chem> getOutputs() {
        return outputs;
    }

    public List<Rxn> getReactions() {
        return reactions;
    }

    public Map<String, Chem> getChemicals() {
        return chemicals;
    }
}
