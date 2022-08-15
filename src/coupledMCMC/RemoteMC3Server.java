package coupledMCMC;

import java.net.*;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.tools.Application;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.Logger;
import beast.base.inference.Runnable;
import beast.base.core.Log;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;

import java.io.*;

@Description("Server for running MC3 remotely: for every chain (hot or cold) one server should be running")
public class RemoteMC3Server extends Runnable {
	final public Input<Integer> portInput = new Input<>("port", "port for this MC3 server to listen to",
			Validate.REQUIRED);

	// initialize socket and input stream
	private Socket socket = null;
	private ServerSocket server = null;
	private DataInputStream in = null;
	private DataOutputStream out = null;

	@Override
	public void initAndValidate() {
	}

	HeatedChain mc3;
	
	@Override
	public void run() {

		// starts server and waits for a connection
		try {
			server = new ServerSocket(portInput.get());
			System.out.println("Server started");

			System.out.println("Waiting for a client at " + server + " ...");

			socket = server.accept();
			System.out.println("Client accepted");

			// takes input from the client socket
			in = new DataInputStream(new BufferedInputStream(socket.getInputStream()));
			out = new DataOutputStream(socket.getOutputStream());

			String line = "";

			// reads message from client until "Over" is sent
			while (!line.equals("Over")) {
				try {
					line = in.readUTF();
					if (line.startsWith("<?xml")) {
						// parse the XML that is sent through
						XMLParser parser = new XMLParser();
						Object o = null;
						try {
							o = parser.parseFragment(line, true);
						} catch (XMLParserException e) {
							e.printStackTrace();
							out.writeUTF("XML parsing failed "+ e.getMessage());
						}
						// make sure it is a HeatedChain object
						if (o instanceof HeatedChain) {
							mc3 = (HeatedChain) o;							
							// remove all screen loggers
							for (int j = mc3.coupledLoggersInput.get().size()-1; j >=0 ; j--){
								if (mc3.coupledLoggersInput.get().get(j).getID().contentEquals("screenlog")){
									mc3.coupledLoggersInput.get().remove(j);
								}
							}
							
							// remove all loggers of heated chains if they are not logged
							for (int j = 0; j < mc3.coupledLoggersInput.get().size(); j++) {
									mc3.coupledLoggersInput.get().get(j).setSuppressLogging(true);
							}
							// initialize each chain individually
							// mc3.setResampleEvery(resampleEvery);
							// mc3.setTemperature(i, getTemperature(i));
							
							// needed to avoid error of putting the working dir twice
							// mc3.setChainNr(i);
							out.writeUTF("mc3 created OK");
						} else {
							out.writeUTF("Expected MCMC wth HeatedChain as main object");
						}
					} else if (line.startsWith("sets")) {
						String stateFileName = line.split(",")[1];
						mc3.setStateFile(stateFileName, restoreFromFile);
						out.writeUTF("state file set to " + mc3.getStateFileName());
					} else if (line.startsWith("setb")) {
						double beta = Double.parseDouble(line.split(",")[1]);
						mc3.setBeta(beta);
						out.writeUTF("beta set to " + beta);
					} else if (line.startsWith("geto")) {
						line = getOperatorParameters();
						out.writeUTF(line);
					} else if (line.startsWith("getu")) {						
						out.writeUTF(mc3.getScaledLogLikelihood(1.0) + "");
					} else if (line.startsWith("seto")) {
						setOperatorParameters(line);
						out.writeUTF("operators updates");
					} else if (line.startsWith("setr")) {
						int resampleEvery = Integer.parseInt(line.split(",")[1]);
						mc3.setResampleEvery(resampleEvery);
						try {
							mc3.run();
						} catch (SAXException | ParserConfigurationException e) {
							e.printStackTrace();
						}
						out.writeUTF("resampleEvery set to " + resampleEvery);
					} else if (line.startsWith("runtill")) {
						int i = Integer.parseInt(line.split(",")[1]);
						try {
							mc3.runTillResample(i);
						} catch (Exception e) {
							out.writeUTF("Failed " + e.getMessage());
						}
						out.writeUTF(mc3.getCurrentLogLikelihood()+"");
					} else if (line.startsWith("seti")) {
						int i = Integer.parseInt(line.split(",")[1]);
						mc3.setChainNr(i);
						// only make the cold chain log
						for (CoupledLogger logger : mc3.coupledLoggersInput.get()) {
							logger.setSuppressLogging(i != 0);
						}
						out.writeUTF("chain nr set to " + i);						
					} else if (line.startsWith("resume")) {
						String resume = line.split("=")[1];
						restoreFromFile = Boolean.parseBoolean(resume);
					}
					
					Log.debug.println(line);
					Log.debug("Message received " + line);
				} catch (IOException e) {
					System.out.println(e);
					break;
				}
			}
			System.out.println("Closing connection");

			// close connection
			socket.close();
			in.close();
		} catch (IOException i) {
			System.out.println(i);
		}
	}

	private String getOperatorParameters() {
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < mc3.operatorsInput.get().size(); i++){
			Operator operator = mc3.operatorsInput.get().get(i);

		    int nrRejected = operator.get_m_nNrRejected();
		    b.append(',').append(nrRejected);
		    int nrAccepted = operator.get_m_nNrAccepted();
		    b.append(',').append(nrAccepted);
		    int nrRejectedForCorrection = operator.get_m_nNrRejectedForCorrection();
		    b.append(',').append(nrRejectedForCorrection);
		    int nrAcceptedForCorrection = operator.get_m_nNrAcceptedForCorrection();
		    b.append(',').append(nrAcceptedForCorrection);
		    
		    double coercableParameterValue = operator.getCoercableParameterValue();
		    b.append(',').append(coercableParameterValue);
			
		}
		
		if (b.toString().split(",").length != mc3.operatorsInput.get().size() * 5 + 1) {
			System.err.println("something is wrong with the length of generated OperatorParameters");
		}
		return b.toString();
	}

	private void setOperatorParameters(String line) {
		String [] strs = line.split(",");
		if (strs.length != mc3.operatorsInput.get().size() * 5 + 1) {
			System.err.println("something is wrong with the length of transmitted OperatorParameters");
		}
		int k = 1;
		for (int i = 0; i < mc3.operatorsInput.get().size(); i++){
			Operator operator = mc3.operatorsInput.get().get(i);

		    int nrRejected = Integer.parseInt(strs[k++]);
		    int nrAccepted = Integer.parseInt(strs[k++]);
		    int nrRejectedForCorrection = Integer.parseInt(strs[k++]);
		    int nrAcceptedForCorrection = Integer.parseInt(strs[k++]);
		    operator.setAcceptedRejected(nrAccepted, nrRejected, nrAcceptedForCorrection, nrRejectedForCorrection);
		    
		    double coercableParameterValue = Double.parseDouble(strs[k++]);
		    operator.setCoercableParameterValue(coercableParameterValue);
		}			
		
	}

	public static void main(String args[]) throws Exception {
		new Application(new RemoteMC3Server(), "Remote MC3 Server", args);
	}
}
