package coupledMCMC;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashSet;
import java.util.Set;

import beast.base.core.Log;

public class RemoteLoggerServer {

	// server is listening on port
	private ServerSocket ss;
	private Thread serverThread;
	Set<String> files;
	boolean restoreFromFile;

	public RemoteLoggerServer(int port, boolean restoreFromFile) throws IOException {
		this.restoreFromFile = restoreFromFile;
		// server is listening on port
		ss = new ServerSocket(port);
		files = new HashSet<>();
		
		serverThread = new Thread() {
			boolean running = true;

			@Override
			public void run() {
				// running infinite loop for getting client requests
				while (running) {
					Socket s = null;

					try {
						// socket object to receive incoming client requests
						s = ss.accept();

						Log.info("A new client is connected : " + s);

						// obtaining input and out streams
						DataInputStream dis = new DataInputStream(s.getInputStream());

						// create a new thread object
						Thread t = new LoggerClientHandler(s, dis);

						// Invoking the start() method
						t.start();

					} catch (Exception e) {
						try {
							s.close();
						} catch (IOException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						e.printStackTrace();
					}
				}
			}
			
			@Override
			public void interrupt() {
				running = false;
				super.interrupt();
			}
		};
			
		serverThread.start();

	}

	public class LoggerClientHandler extends Thread {
		final DataInputStream dis;
		final Socket s;
		PrintStream out = null;

		// Constructor
		public LoggerClientHandler(Socket s, DataInputStream dis) {
			this.s = s;
			this.dis = dis;
		}

		@Override
		public void run() {

			// first message indicates the file name
			try {
				String file = dis.readUTF();
				if (file.equals("null")) {
					out = System.out;
				} else {
					// TODO: deal with resuming
					if (new File(file).exists() || files.contains(file)) {
						Log.warning("Appending file " + file);
						out = new PrintStream(new FileOutputStream(file, true));
					} else {
						Log.warning("Writing file " + file);
						out = new PrintStream(file);
						files.add(file);
					}
				}
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			String received;
			while (true) {
				try {

					// receive the answer from client
					received = dis.readUTF();
					if (received.equals("close")) {
						break;
					}
					if (received.length() > 0) {
						out.print(received);
						out.flush();
					}
				} catch (EOFException e) {
					break;
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			try {
				// closing resources
				this.dis.close();
				// close connections
				serverThread.interrupt();
				ss.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
