package beast.coupledMCMC;

import java.io.PrintStream;

import beast.core.Logger;

public class CoupledLogger extends Logger {

    public void setM_out(PrintStream m_out_alt) {
    	m_out = m_out_alt;  
    }

}
