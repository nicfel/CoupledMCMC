package coupledMCMC;


import java.io.FileWriter;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.XMLFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.inference.Runnable;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;

@Description("Convert MCMC analysis to a coupled MCMC analysis")
public class MCMC2CoupledMCMC extends Runnable {
	
	public Input<XMLFile> model1Input = new Input<>("xml",
			"file name of BEAST XML file containing the model for which to create a coupled MCMC XML file for",
			new XMLFile("examples/normalTest-1XXX.xml"), Validate.REQUIRED);
	public Input<OutFile> outputInput = new Input<>("output", "where to save the file", new OutFile("beast.xml"));
	
	
	public Input<Integer> nrOfChainsInput = new Input<Integer>("chains", " number of chains to run in parallel (default 2)", 2);
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 100);
	public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written","");
	
	// input of the difference between temperature scalers
	public Input<Double> deltaTemperatureInput = new Input<>("deltaTemperature","temperature difference between the i-th and the i-th+1 chain", 0.1);
	public Input<Double> maxTemperatureInput = new Input<>("maxTemperature","maximal temperature of the hottest chain");	
//	public Input<Boolean> logHeatedChainsInput = new Input<>("logHeatedChains","if true, log files for heated chains are also printed", false);
	
	// Input on whether the temperature between chains should be optimized
	public Input<Boolean> optimiseInput = new Input<>("optimise","if true, the temperature is automatically optimised to reach a target acceptance probability", true);
	
	public Input<Integer> optimiseDelayInput = new Input<>("optimiseDelay","after this many epochs/swaps, the temperature will be optimized (if optimising is set to true)", 100);
	public Input<Double> targetAcceptanceProbabilityInput = new Input<>("target", "target acceptance probability of swaps", 0.234);
	
//	public Input<Boolean> preScheduleInput = new Input<>("preSchedule","if true, how long chains are run for is scheduled at the beginning", true);

	
	@Override
	public void initAndValidate() {
	}

	
	@Override
	public void run() {
		XMLParser parser = new XMLParser();
		MCMC mcmc;
		try {
			mcmc = (MCMC) parser.parseFile(model1Input.get());
			
			CoupledMCMC mc3 = new CoupledMCMC();
			// nested sampling options
			mc3.nrOfChainsInput.setValue(nrOfChainsInput.get(), mc3);
			mc3.resampleEveryInput.setValue(resampleEveryInput.get(), mc3);
			mc3.tempDirInput.setValue(tempDirInput.get(), mc3);
			
			mc3.deltaTemperatureInput.setValue(deltaTemperatureInput.get(), mc3);
			mc3.maxTemperatureInput.setValue(maxTemperatureInput.get(), mc3);
			mc3.optimiseInput.setValue(optimiseInput.get(), mc3);
			
			mc3.logHeatedChainsInput.setValue(false, mc3);
			mc3.preScheduleInput.setValue(true, mc3);

			// mcmc options
			mc3.posteriorInput.setValue(mcmc.posteriorInput.get(), mc3);
			mc3.startStateInput.setValue(mcmc.startStateInput.get(), mc3);
			mc3.operatorsInput.setValue(mcmc.operatorsInput.get(), mc3);
			mc3.chainLengthInput.setValue(mcmc.chainLengthInput.get(), mc3);
			mc3.initialisersInput.setValue(mcmc.initialisersInput.get(), mc3);			
			mc3.loggersInput.setValue(mcmc.loggersInput.get(), mc3);
			mc3.operatorScheduleInput.setValue(mcmc.operatorScheduleInput.get(), mc3);

		
	        Log.warning("Writing to file " + outputInput.get().getPath());
			XMLProducer producer = new XMLProducer();
			String xml = producer.toXML(mc3);
	        FileWriter outfile = new FileWriter(outputInput.get());
	        outfile.write(xml);
	        outfile.close();
	        
			
		} catch (SAXException | IOException | ParserConfigurationException | XMLParserException e) {
			System.out.println(e);
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
        Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new MCMC2CoupledMCMC(), "Convert MCMC to Coupled MCMC", args);
	}

}
