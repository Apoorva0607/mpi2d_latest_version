setpref('Internet', 'E_mail', 'devavrat0210@gmail.com');
setpref('Internet', 'SMTP_Username', 'devavrat0210@gmail.com');
setpref('Internet', 'SMTP_Password', 'devil@ndroid0210');
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port', '465');

sendmail('devavrat.likhite@utah.edu','Done!!!!',...
'Done');