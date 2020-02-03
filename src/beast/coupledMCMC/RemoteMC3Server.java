package beast.coupledMCMC;

import java.net.*;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.Runnable;
import beast.util.XMLParser;
import beast.util.XMLParserException;

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
							// String[] splittedFileName = stateFileName.split("/");
							
							// mc3.setStateFile(
							//		splittedFileName[splittedFileName.length-1].replace(".state", "." + i + "state"), restoreFromFile);
							// mc3.setChainNr(i);
							out.writeUTF("mc3 created OK");
						} else {
							out.writeUTF("Expected MCMC wth HeatedChain as main object");
						}
					} else if (line.startsWith("setb")) {
						double beta = Double.parseDouble(line.split(",")[1]);
						mc3.setBeta(beta);
						out.writeUTF("beta set to " + beta);
					} else if (line.startsWith("geto")) {
						line = getOperatorParameters();
						out.writeUTF(line);
					} else if (line.startsWith("getu")) {						
						out.writeUTF(mc3.getUnscaledCurrentLogLikelihood() + "");
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
						for (CoupledLogger logger : mc3.loggers) {
							logger.setSuppressLogging(i == 0);
						}
						out.writeUTF("chain nr set to " + i);						
					}
					
					System.out.println(line);
					out.writeUTF("Message received " + line);
				} catch (IOException i) {
					System.out.println(i);
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
			Operator operator1 = mc3.operatorsInput.get().get(i);

		    int m_nNrRejected = operator1.get_m_nNrRejected();
		    b.append(',').append(m_nNrRejected);
		    int m_nNrAccepted = operator1.get_m_nNrAccepted();
		    b.append(',').append(m_nNrAccepted);
		    int m_nNrRejectedForCorrection = operator1.get_m_nNrRejectedForCorrection();
		    b.append(',').append(m_nNrRejectedForCorrection);
		    int m_nNrAcceptedForCorrection = operator1.get_m_nNrAcceptedForCorrection();
		    b.append(',').append(m_nNrAcceptedForCorrection);
		    
		    double coercableParameterValue = operator1.getCoercableParameterValue();
		    b.append(',').append(coercableParameterValue);
			
		}			
		return b.toString();
	}

	private void setOperatorParameters(String line) {
		String [] strs = line.split(",");
		int k = 1;
		for (int i = 0; i < mc3.operatorsInput.get().size(); i++){
			Operator operator1 = mc3.operatorsInput.get().get(i);

		    int m_nNrRejected = Integer.parseInt(strs[k++]);
		    int m_nNrAccepted = Integer.parseInt(strs[k++]);
		    int m_nNrRejectedForCorrection = Integer.parseInt(strs[k++]);
		    int m_nNrAcceptedForCorrection = Integer.parseInt(strs[k++]);
		    
		    operator1.setAcceptedRejected(m_nNrAccepted, m_nNrRejected, m_nNrAcceptedForCorrection, m_nNrRejectedForCorrection);
		    
		    
		    double coercableParameterValue = Double.parseDouble(strs[k++]);
		    operator1.setCoercableParameterValue(coercableParameterValue);
		}			
		
	}

	public static void main(String args[]) throws Exception {
		new Application(new RemoteMC3Server(), "Remote MC3 Server", args);
	}
}
