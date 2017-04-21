package japsafuture.stream;

import static java.lang.Thread.currentThread;
import static java.util.concurrent.Executors.newSingleThreadExecutor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Flow.Publisher;
import java.util.concurrent.Flow.Subscriber;
import java.util.concurrent.Flow.Subscription;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class MyPublisherSimple implements Publisher<Integer> {

    private static final String LOG_MESSAGE_FORMAT = "Publisher >> [%s] %s%n";
    private List<MySubscription> subscriptions = Collections.synchronizedList(new ArrayList<MySubscription>());

    @Override
    public void subscribe(Subscriber<? super Integer> subscriber) {
        MySubscription subscription = new MySubscription(subscriber);

        subscriptions.add(subscription);
        subscriber.onSubscribe(subscription);
    }

    private class MySubscription implements Subscription {

        private Subscriber<? super Integer> subscriber;
        private final AtomicInteger value;
        private AtomicBoolean isCanceled;

        public MySubscription(Subscriber<? super Integer> subscriber) {
            this.subscriber = subscriber;

            value = new AtomicInteger();
            isCanceled = new AtomicBoolean(false);
        }

        @Override
        public void request(long n) {
            if (isCanceled.get())
                return;

            publishItems(n);
        }

        @Override
        public void cancel() {
            isCanceled.set(true);

            synchronized (subscriptions) {
                subscriptions.remove(this);
                if (subscriptions.size() == 0)
                    shutdown();
            }
        }

        private void publishItems(long n) {
            for (int i = 0; i < n; i++) {
                int v = value.incrementAndGet();
                   log("publish item: [" + v + "] ...");
                subscriber.onNext(v);

            }
        }

        private void shutdown() {
            log("Shut down executor...");
        }

    }

    private void log(String message, Object... args) {
        String fullMessage = String.format(LOG_MESSAGE_FORMAT, currentThread().getName(), message);

        System.out.printf(fullMessage, args);
    }
}