package beast.coupledMCMC;

import java.io.IOException;
import java.io.PrintStream;

import beast.core.Logger;

public class CoupledLogger extends Logger {
	
	public CoupledLogger(Logger other) {
		loggersInput.setValue(other.loggersInput.get(), this);
		fileNameInput.setValue(other.fileNameInput.get(), this);
		everyInput.setValue(other.everyInput.get(), this);
		modelInput.setValue(other.modelInput.get(), this);
		modeInput.setValue(other.modeInput.get(), this);
		sortModeInput.setValue(other.sortModeInput.get(), this);
		sanitiseHeadersInput.setValue(other.sanitiseHeadersInput.get(), this);
	}

	public CoupledLogger() {
	}
	
	boolean suppressLogging = false;
	
    public void setM_out(PrintStream m_out_alt) {
    	m_out = m_out_alt;  
    }
	
    public void setSuppressLogging(boolean suppressLogging) {
    	this.suppressLogging = suppressLogging;  
    }
    
    public boolean getSuppressLogging() {
    	return suppressLogging;  
    }

    @Override
    public void init() throws IOException{
    	if (!suppressLogging){
    		super.init();
    	}
    }

    
    
//    @Override
//	protected boolean openLogFile() throws IOException {
//    	if (!suppressLogging)
//    		super.openLogFile();
//    	
//    	return false;
//	}
    
    @Override
    public void log(long sampleNr){
    	if (!suppressLogging)
    		super.log(sampleNr);
    }
    
    @Override
    public void close() {
    	if (!suppressLogging)
    		super.close();
    }
}
