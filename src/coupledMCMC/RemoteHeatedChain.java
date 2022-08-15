package coupledMCMC;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.Socket;
import java.net.UnknownHostException;

import beast.base.core.Description;
import beast.base.core.Log;

@Description("Local representation of heated chain running on a remote machine")
public class RemoteHeatedChain {
	// initialize socket and input output streams
	private Socket socket = null;
	private DataInputStream sin = null;
	private DataOutputStream out = null;

	private double beta;
	private int chainNr = 0;
	private double oldLogLikelihood = 0;

	public RemoteHeatedChain(String xml, String host, int port, boolean restoreFromFile) {

		// establish a connection
		try {
			socket = new Socket(host, port);
			System.out.println("Connected");

			// sends output to the socket
			out = new DataOutputStream(socket.getOutputStream());

			// reads responses from the socket
			sin = new DataInputStream(socket.getInputStream());
			out.writeUTF("resume=" + restoreFromFile);
			out.writeUTF(xml);
			String response = sin.readUTF();
			System.err.println(response);
		} catch (UnknownHostException u) {
			System.out.println(u);
		} catch (IOException i) {
			System.out.println(i);
		}

	}

	protected double getCurrentLogLikelihood() {
		return oldLogLikelihood * beta;
	}

	protected double getUnscaledCurrentLogLikelihood() {
		try {
			out.writeUTF("getu");
			String response = sin.readUTF();
			log(response);
			return Double.parseDouble(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return Double.NaN;
	}

	protected void setResampleEvery(int resampleEvery) {
		try {
			out.writeUTF("setr," + resampleEvery);
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
		// this.resampleEvery = resampleEvery;
	}

	// // RRB: what is the purpose of the first argument?
	protected void setTemperature(int i, double temperature) {
		try {
			out.writeUTF("setb," + 1 / (1 + temperature));
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.beta = 1 / (1 + temperature);
	}

	protected double getBeta() {
		return beta;
	}

	protected void setBeta(double beta) {
		try {
			out.writeUTF("setb," + beta);
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.beta = beta;
	}

	protected int getChainNr() {
		return chainNr;
	}

	private void log(String s) {
		// Log.debug(s);
	}
	
	protected void setChainNr(int chainNr) {
		this.chainNr = chainNr;
		try {
			out.writeUTF("seti," + chainNr);
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	protected void setStateFile(String stateFileName, boolean restoreFromFile) {
		try {
			out.writeUTF("sets," + stateFileName + "," + restoreFromFile);
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// run MCMC inner loop for resampleEvery nr of samples
	protected long runTillResample(long runUntil) throws Exception {
		try {
			out.writeUTF("runtill," + runUntil);
			String response = sin.readUTF();
			log(response);
			oldLogLikelihood = Double.parseDouble(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return System.currentTimeMillis();
	}

	public String getOperatorTuningInfo() {
		try {
			out.writeUTF("geto");
			String response = sin.readUTF();
			log(response);
			return response;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	public void setOperatorTuningInfo(String tuning) {
		try {
			out.writeUTF("seto" + tuning);
			String response = sin.readUTF();
			log(response);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
