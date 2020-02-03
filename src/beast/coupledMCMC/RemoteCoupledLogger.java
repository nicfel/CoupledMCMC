package beast.coupledMCMC;

import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.InetAddress;
import java.net.Socket;
import java.net.UnknownHostException;
import java.nio.charset.StandardCharsets;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;

@Description("Logger for CoupledMCMC that logs via a RemoteLoggerService")
class RemoteCoupledLogger extends CoupledLogger {
	final public Input<String> hostInput = new Input<>("host", "URL of MC3 logger service", Validate.REQUIRED);
	final public Input<Integer> portInput = new Input<>("port", "port of MC3 logger service", Validate.REQUIRED);

	private Socket socket = null;
	private DataOutputStream out = null;

	@Override
	public void initAndValidate() {

		// establish a connection
		try {
			String address = InetAddress.getByName(hostInput.get()).getHostAddress();
			socket = new Socket(address, portInput.get());
			Log.info("Connected to logger service at " + address + ":" + portInput.get());
		
			// sends output to the socket
			out = new DataOutputStream(socket.getOutputStream());

		} catch (UnknownHostException u) {
			System.out.println(u);
		} catch (IOException i) {
			System.out.println(i);
		}

		super.initAndValidate();
	}

	@Override
	public void init() throws IOException {
		if (fileNameInput.get() != null) {
			out.writeUTF(fileNameInput.get());
		} else {
			out.writeUTF("null");
		}

		if (!suppressLogging) {
			final ByteArrayOutputStream baos = new ByteArrayOutputStream();
			final String utf8 = StandardCharsets.UTF_8.name();
			try {
				m_out = new PrintStream(baos, true, utf8);
				super.init();
				String data = baos.toString(utf8);
				out.writeUTF(data);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Override
	public void log(long sampleNr) {
		if (!suppressLogging) {
			final ByteArrayOutputStream baos = new ByteArrayOutputStream();
			final String utf8 = StandardCharsets.UTF_8.name();
			try {
				m_out = new PrintStream(baos, true, utf8);
				super.log(sampleNr);
				String data = baos.toString(utf8);
				out.writeUTF(data);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Override
    public void close() {
    	if (!suppressLogging) {
    	    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    	    final String utf8 = StandardCharsets.UTF_8.name();    	    
    	    try {
				m_out = new PrintStream(baos, true, utf8);
	    		super.close();
	    		String data = baos.toString(utf8);
				out.writeUTF(data);
				out.writeUTF("close");
			} catch (IOException e) {
				e.printStackTrace();
			}
    	}
    }

}
